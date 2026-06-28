# Copyright 2024 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Client API."""

from __future__ import annotations

import datetime
import gzip
import json
import os
import re
from collections.abc import AsyncIterator, Awaitable
from contextlib import asynccontextmanager
from dataclasses import asdict, dataclass
from typing import Any, cast
from uuid import uuid4

import psycopg
from fastapi import (
    APIRouter,
    BackgroundTasks,
    Depends,
    HTTPException,
    Query,
    Response,
    status,
)
from ord_schema.logging import get_logger
from ord_schema.orm.database import get_connection_string
from ord_schema.proto import dataset_pb2
from psycopg import AsyncCursor
from psycopg.rows import dict_row
from pydantic import BaseModel
from rdkit import Chem
from redis.asyncio import Redis

from ord_interface.api.queries import (
    DatasetIdQuery,
    DoiQuery,
    QueryResult,
    ReactionComponentQuery,
    ReactionConversionQuery,
    ReactionIdQuery,
    ReactionSmartsQuery,
    ReactionYieldQuery,
    StatsResult,
    fetch_dataset_most_used_smiles_for_inputs,
    fetch_dataset_most_used_smiles_for_products,
    fetch_reactions,
    run_queries,
)

logger = get_logger(__name__)
router = APIRouter(tags=["client"])

BOND_LENGTH = 20
MAX_RESULTS = 1000


@asynccontextmanager
async def get_cursor() -> AsyncIterator[AsyncCursor[dict[str, Any]]]:
    """Returns a psycopg cursor."""
    dsn = os.getenv("ORD_INTERFACE_POSTGRES")
    if dsn is None:
        dsn = re.sub(
            r"postgresql\+psycopg",
            "postgresql",
            get_connection_string(
                database=os.environ["POSTGRES_DATABASE"],
                username=os.environ["POSTGRES_USER"],
                password=os.environ["POSTGRES_PASSWORD"],
                host=os.environ["POSTGRES_HOST"],
            ),
        )
    async with await psycopg.AsyncConnection[dict[str, Any]].connect(
        dsn, row_factory=dict_row, options="-c search_path=public,ord"
    ) as connection:
        await connection.set_read_only(True)
        async with connection.cursor() as cursor:
            yield cursor


@asynccontextmanager
async def get_redis() -> AsyncIterator[Redis]:
    """Returns a Redis client instance."""
    host = os.environ.get("REDIS_HOST", "localhost")
    port = int(os.environ.get("REDIS_PORT", "6379"))
    ssl = os.environ.get("REDIS_SSL", "0") == "1"
    async with Redis(host=host, port=port, ssl=ssl) as client:
        # redis.asyncio stubs return Awaitable[bool] | bool from ping(); the runtime is always awaitable.
        if not await cast(Awaitable[bool], client.ping()):
            raise RuntimeError(
                f"Failed to connect to Redis server {host}:{port} ({ssl=})"
            )
        logger.debug(f"Connected to Redis server {host}:{port} ({ssl=})")
        yield client


@dataclass
class QueryParams:
    """Query parameters."""

    # NOTE(skearnes): BaseModel does not work here; see https://github.com/fastapi/fastapi/discussions/10556.

    dataset_id: list[str] | None = Query(None)
    reaction_id: list[str] = Query(None)
    reaction_smarts: str | None = None
    min_conversion: float | None = None
    max_conversion: float | None = None
    min_yield: float | None = None
    max_yield: float | None = None
    doi: list[str] | None = Query(None)
    # Each value is a JSON-encoded ComponentSpec; see ComponentSpec.
    component: list[str] | None = Query(None)
    use_stereochemistry: bool | None = None
    similarity: float | None = None
    limit: int | None = None


class ComponentSpec(BaseModel):
    """A single component predicate, JSON-encoded in each ``component`` query value.

    Replaces a legacy ``"pattern;target;mode"`` string whose ``;`` delimiter collided
    with SMARTS patterns (which use ``;`` as a low-precedence AND). ``target`` and
    ``mode`` are matched case-insensitively against the ReactionComponentQuery enums.
    """

    pattern: str
    target: str
    mode: str

    @classmethod
    def parse(cls, spec: str) -> ComponentSpec:
        """Parses one ``component`` value, accepting JSON or the legacy format.

        New values are JSON objects; the legacy ``"pattern;target;mode"`` form is still
        accepted so previously shared search URLs keep working. The legacy parse splits
        from the right because a SMARTS ``pattern`` may itself contain ``;``.

        Args:
            spec: A single ``component`` query-parameter value.

        Returns:
            The parsed ComponentSpec.

        Raises:
            ValueError: If ``spec`` is neither valid JSON nor a 3-field legacy string.
        """
        spec = spec.strip()
        if spec.startswith("{"):
            return cls.model_validate_json(spec)
        pattern, target, mode = spec.rsplit(";", 2)
        return cls(pattern=pattern, target=target, mode=mode)


async def run_query(
    params: QueryParams, return_ids: bool
) -> list[QueryResult] | list[str]:
    """Runs a query and returns a list of matched reactions."""
    queries = []
    if params.dataset_id and isinstance(params.dataset_id, list):
        queries.append(DatasetIdQuery(params.dataset_id))
    if params.reaction_id and isinstance(params.reaction_id, list):
        queries.append(ReactionIdQuery(params.reaction_id))
    if params.reaction_smarts:
        queries.append(ReactionSmartsQuery(params.reaction_smarts))
    if params.min_conversion is not None or params.max_conversion is not None:
        queries.append(
            ReactionConversionQuery(params.min_conversion, params.max_conversion)
        )
    if params.min_yield is not None or params.max_yield is not None:
        queries.append(ReactionYieldQuery(params.min_yield, params.max_yield))
    if params.doi and isinstance(params.doi, list):
        queries.append(DoiQuery(params.doi))
    if params.component and isinstance(params.component, list):
        kwargs = {}
        if params.use_stereochemistry is not None:
            kwargs["use_chirality"] = params.use_stereochemistry
        if params.similarity is not None:
            kwargs["similarity_threshold"] = params.similarity
        for spec in params.component:
            # Cover the whole predicate: a bad target/mode is a KeyError on the enum
            # lookup and an unparseable pattern is a ValueError from the query -- both
            # are client errors, not 500s.
            try:
                component = ComponentSpec.parse(spec)
                queries.append(
                    ReactionComponentQuery(
                        component.pattern,
                        ReactionComponentQuery.Target[component.target.upper()],
                        ReactionComponentQuery.MatchMode[component.mode.upper()],
                        **kwargs,
                    )
                )
            except (KeyError, ValueError) as error:
                raise HTTPException(
                    status_code=400, detail=f"Invalid component spec: {spec!r}"
                ) from error
    if not queries:
        raise ValueError("No query parameters were specified.")
    limit = MAX_RESULTS
    if params.limit:
        limit = min(params.limit, MAX_RESULTS)
    async with get_cursor() as cursor:
        reaction_ids = await run_queries(cursor, queries, limit=limit)
        if return_ids:
            return reaction_ids
        return await fetch_reactions(cursor, reaction_ids)


@router.get("/query")
async def query(params: QueryParams = Depends()) -> list[QueryResult]:
    """Runs a query."""
    result = await run_query(params, return_ids=False)
    return cast(list[QueryResult], result)  # Type hint.


@router.get("/reaction")
async def get_reaction(reaction_id: str) -> QueryResult:
    """Fetches a Reaction by ID."""
    async with get_cursor() as cursor:
        results = await fetch_reactions(cursor, [reaction_id])
    return results[0]


class ReactionIdList(BaseModel):
    """Reaction ID input."""

    reaction_ids: list[str]


@router.post("/reactions")
async def get_reactions(inputs: ReactionIdList) -> list[QueryResult]:
    """Fetches a list of Reactions by ID."""
    async with get_cursor() as cursor:
        return await fetch_reactions(cursor, inputs.reaction_ids)


class DatasetInfo(BaseModel):
    """Dataset info."""

    dataset_id: str
    name: str
    description: str | None
    num_reactions: int
    submitted_at: datetime.date | None


@router.get("/datasets")
async def get_datasets() -> list[DatasetInfo]:
    """Returns info about the current datasets, most recently submitted first."""
    async with get_cursor() as cursor:
        await cursor.execute(
            "SELECT dataset_id, name, description, num_reactions, submitted_at FROM dataset "
            "ORDER BY submitted_at DESC NULLS LAST, dataset_id"
        )
        return [DatasetInfo(**row) async for row in cursor]


@router.get("/dataset")
async def get_dataset(dataset_id: str) -> DatasetInfo:
    """Returns info about a single dataset."""
    async with get_cursor() as cursor:
        await cursor.execute(
            "SELECT dataset_id, name, description, num_reactions, submitted_at FROM dataset WHERE dataset_id = %s",
            (dataset_id,),
        )
        row = await cursor.fetchone()
        if row is None:
            raise HTTPException(
                status_code=404, detail=f"dataset not found: {dataset_id}"
            )
        return DatasetInfo(**row)


@router.get("/molfile")
async def get_molfile(smiles: str) -> str:
    """Returns a molblock for the given SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(status_code=400, detail=f"invalid SMILES: {smiles}")
    return Chem.MolToMolBlock(mol)


@router.post("/download_search_results")
async def get_search_results(inputs: ReactionIdList):
    """Downloads search results as a Dataset proto."""
    results = await get_reactions(inputs)
    dataset = dataset_pb2.Dataset(
        name="ORD Search Results", reactions=[result.reaction for result in results]
    )
    return Response(
        gzip.compress(dataset.SerializeToString()), media_type="application/gzip"
    )


async def run_task(task_id: str, params: QueryParams) -> bool:
    """Wraps run_query() as a background task."""
    # NOTE(skearnes): Use reaction IDs to avoid stuffing full protos into the result database.
    result = await run_query(params, return_ids=True)
    logger.debug(f"Finished task {task_id}")
    async with get_redis() as client:
        return await client.set(f"result:{task_id}", json.dumps(result), ex=60 * 60)


@router.get("/submit_query")
async def submit_query(
    background_tasks: BackgroundTasks, params: QueryParams = Depends()
) -> str:
    """Submits a query as a background task."""
    task_id = str(uuid4())
    async with get_redis() as client:
        await client.set(f"query:{task_id}", json.dumps(asdict(params)), ex=60 * 60)
    background_tasks.add_task(run_task, task_id=task_id, params=params)
    logger.debug(f"Created task {task_id}")
    return task_id


@router.get("/fetch_query_result")
async def fetch_query_result(task_id: str):
    """Checks the query status, returning the results if the query is complete."""
    async with get_redis() as client:
        if not await client.exists(f"query:{task_id}"):
            return Response(
                f"Task {task_id} does not exist", status_code=status.HTTP_404_NOT_FOUND
            )
        result = await client.get(f"result:{task_id}")
    if result is None:
        return Response(
            f"Task {task_id} is pending", status_code=status.HTTP_202_ACCEPTED
        )
    async with get_cursor() as cursor:
        return await fetch_reactions(cursor, json.loads(result))


@router.get("/input_stats")
async def get_input_stats(dataset_id: str, limit: int = 30) -> list[StatsResult]:
    """Return the most frequently used input SMILES for a dataset."""
    async with get_cursor() as cursor:
        results = await fetch_dataset_most_used_smiles_for_inputs(
            cursor, dataset_id, limit=limit
        )
    return results


@router.get("/product_stats")
async def get_product_stats(dataset_id: str, limit: int = 30) -> list[StatsResult]:
    """Return the most frequently used product SMILES for a dataset."""
    async with get_cursor() as cursor:
        results = await fetch_dataset_most_used_smiles_for_products(
            cursor, dataset_id, limit=limit
        )
    return results
