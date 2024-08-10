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

import gzip
import os
import re
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from typing import Iterator

import psycopg
from celery import shared_task
from celery.result import AsyncResult
from celery.states import SUCCESS, UNREADY_STATES
from fastapi import APIRouter, Depends, Query, Response, status
from ord_schema.orm.database import get_connection_string
from ord_schema.proto import dataset_pb2
from psycopg import Cursor
from psycopg.rows import dict_row
from pydantic import BaseModel
from rdkit import Chem

from ord_interface.api.queries import (
    DatasetIdQuery,
    DoiQuery,
    QueryResult,
    ReactionComponentQuery,
    ReactionConversionQuery,
    ReactionIdQuery,
    ReactionSmartsQuery,
    ReactionYieldQuery,
    run_queries,
)

router = APIRouter(prefix="/api", tags=["client"])

BOND_LENGTH = 20
MAX_RESULTS = 1000


@contextmanager
def get_cursor() -> Iterator[Cursor]:
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
    with psycopg.connect(  # pylint: disable=not-context-manager
        dsn, row_factory=dict_row, options="-c search_path=public,ord"
    ) as connection:
        connection.set_read_only(True)
        with connection.cursor() as cursor:
            yield cursor


@dataclass
class QueryParams:  # pylint: disable=too-many-instance-attributes
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
    component: list[str] | None = Query(None)
    use_stereochemistry: bool | None = None
    similarity: float | None = None
    limit: int | None = None


@router.get("/query")
def query(params: QueryParams = Depends()) -> list[QueryResult]:
    """Runs a query and returns a list of matched reactions."""
    # pylint: disable=too-many-arguments,too-many-branches,too-many-locals
    queries = []
    if params.dataset_id and isinstance(params.dataset_id, list):
        queries.append(DatasetIdQuery(params.dataset_id))
    if params.reaction_id and isinstance(params.reaction_id, list):
        queries.append(ReactionIdQuery(params.reaction_id))
    if params.reaction_smarts:
        queries.append(ReactionSmartsQuery(params.reaction_smarts))
    if params.min_conversion is not None or params.max_conversion is not None:
        queries.append(ReactionConversionQuery(params.min_conversion, params.max_conversion))
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
            pattern, target_name, mode_name = spec.split(";")
            queries.append(
                ReactionComponentQuery(
                    pattern,
                    ReactionComponentQuery.Target[target_name.upper()],
                    ReactionComponentQuery.MatchMode[mode_name.upper()],
                    **kwargs,
                )
            )
    if not queries:
        raise ValueError("No query parameters were specified.")
    limit = MAX_RESULTS
    if params.limit:
        limit = min(params.limit, MAX_RESULTS)
    with get_cursor() as cursor:
        return run_queries(cursor, queries, limit=limit)


@router.get("/reaction")
async def get_reaction(reaction_id: str) -> QueryResult:
    """Fetches a Reaction by ID."""
    with get_cursor() as cursor:
        results = run_queries(cursor, ReactionIdQuery([reaction_id]))
    return results[0]


class ReactionIdList(BaseModel):
    """Reaction ID input."""

    reaction_ids: list[str]


@router.post("/reactions")
async def get_reactions(inputs: ReactionIdList) -> list[QueryResult]:
    """Fetches a list of Reactions by ID."""
    with get_cursor() as cursor:
        return run_queries(cursor, ReactionIdQuery(inputs.reaction_ids))


class DatasetInfo(BaseModel):
    """Dataset info."""

    dataset_id: str
    name: str
    description: str | None
    num_reactions: int


@router.get("/datasets")
async def get_datasets() -> list[DatasetInfo]:
    """Returns info about the current datasets."""
    with get_cursor() as cursor:
        cursor.execute("SELECT dataset_id, name, description, num_reactions FROM dataset")
        return [DatasetInfo(**row) for row in cursor]


@router.get("/molfile")
async def get_molfile(smiles: str) -> str:
    """Returns a molblock for the given SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(smiles)
    return Chem.MolToMolBlock(mol)


@router.post("/download_search_results")
async def get_search_results(inputs: ReactionIdList):
    """Downloads search results as a Dataset proto."""
    results = await get_reactions(inputs)
    dataset = dataset_pb2.Dataset(name="ORD Search Results", reactions=[result.reaction for result in results])
    return Response(gzip.compress(dataset.SerializeToString()), media_type="application/gzip")


class QueryResults(BaseModel):
    """Container for a query results. Only used for celery task result serialization."""

    results: list[QueryResult]


@shared_task(track_started=True)
def run_task(config: dict) -> dict:
    """Wraps query() for celery."""
    results = query(QueryParams(**config))
    return QueryResults(results=results).model_dump()


@router.get("/submit_task")
def submit_task(params: QueryParams = Depends()) -> str:
    """Submits a query as a celery task."""
    task = run_task.delay(asdict(params))
    return task.id


@router.get("/fetch_task")
def fetch_task(task_id: str):
    """Checks the query status, returning the results if the query is complete."""
    task = AsyncResult(task_id)
    state = task.state
    if state in UNREADY_STATES:
        return Response(state, status_code=status.HTTP_102_PROCESSING, media_type="text/plain")
    if state == SUCCESS:
        results = QueryResults(**task.get())
        return results.results
    return Response(state, status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, media_type="text/plain")
