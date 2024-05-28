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

import os
from base64 import b64encode
from typing import Annotated

from fastapi import APIRouter, Query
from pydantic import BaseModel

from ord_interface.client.query import (
    DatasetIdQuery,
    DoiQuery,
    OrdPostgres,
    QueryException,
    ReactionComponentPredicate,
    ReactionComponentQuery,
    ReactionConversionQuery,
    ReactionIdQuery,
    ReactionQuery,
    ReactionSmartsQuery,
    ReactionYieldQuery,
    Result,
)

router = APIRouter(prefix="/api", tags=["client"])

POSTGRES_HOST = os.getenv("POSTGRES_HOST", "localhost")
POSTGRES_PORT = os.getenv("POSTGRES_PORT", "5432")
POSTGRES_USER = os.getenv("POSTGRES_USER", "ord-postgres")
POSTGRES_PASSWORD = os.getenv("POSTGRES_PASSWORD", "ord-postgres")
POSTGRES_DATABASE = os.getenv("POSTGRES_DATABASE", "ord")

BOND_LENGTH = 20
MAX_RESULTS = 1000


def connect():
    return OrdPostgres(
        dbname=POSTGRES_DATABASE,
        user=POSTGRES_USER,
        password=POSTGRES_PASSWORD,
        host=POSTGRES_HOST,
        port=int(POSTGRES_PORT),
    )


def run_query(commands: list[ReactionQuery], limit: int | None) -> list[Result]:
    """Runs a query and returns the matched reactions."""
    if not commands:
        results = []
    elif len(commands) == 1:
        results = connect().run_query(commands[0], limit=limit)
    else:
        # Perform each query without limits and take the intersection of the matched reactions.
        connection = connect()
        reactions = {}
        intersection = None
        for command in commands:
            this_results = {result.reaction_id: result for result in connection.run_query(command)}
            reactions |= this_results
            if intersection is None:
                intersection = set(this_results.keys())
            else:
                intersection &= this_results.keys()
        results = [reactions[reaction_id] for reaction_id in intersection]
        if limit:
            results = results[:limit]
    return results


class QueryResult(BaseModel):
    """Query result."""

    dataset_id: str
    reaction_id: str
    reaction: str | None = None  # Serialized Reaction protocol buffer (base64).

    @classmethod
    def from_result(cls, result: Result) -> QueryResult:
        """Creates a QueryResult from a Result."""
        return cls(
            dataset_id=result.dataset_id,
            reaction_id=result.reaction_id,
            reaction=b64encode(result.proto).decode() if result.reaction else None,
        )

    @classmethod
    def from_results(cls, results: list[Result]) -> list[QueryResult]:
        return [cls.from_result(result) for result in results]


@router.get("/query")
def query(
    dataset_id: Annotated[list[str] | None, Query()] = None,
    reaction_id: Annotated[list[str] | None, Query()] = None,
    reaction_smarts: str | None = None,
    min_conversion: float | None = None,
    max_conversion: float | None = None,
    min_yield: float | None = None,
    max_yield: float | None = None,
    doi: Annotated[list[str] | None, Query()] = None,
    component: Annotated[list[str] | None, Query()] = None,
    use_stereochemistry: bool | None = None,
    similarity: float | None = None,
    limit: int | None = None,
) -> list[QueryResult]:
    """Returns a serialized Dataset proto containing matched reactions."""
    queries = []
    if dataset_id:
        queries.append(DatasetIdQuery(dataset_id))
    if reaction_id:
        queries.append(ReactionIdQuery(reaction_id))
    if reaction_smarts:
        queries.append(ReactionSmartsQuery(reaction_smarts))
    if min_conversion is not None or max_conversion is not None:
        queries.append(ReactionConversionQuery(min_conversion, max_conversion))
    if min_yield is not None or max_yield is not None:
        queries.append(ReactionYieldQuery(min_yield, max_yield))
    if doi:
        queries.append(DoiQuery(doi))
    if component:
        predicates = []
        for spec in component:
            pattern, target_name, mode_name = spec.split(";")
            target = ReactionComponentPredicate.Target.from_name(target_name)
            mode = ReactionComponentPredicate.MatchMode.from_name(mode_name)
            predicates.append(ReactionComponentPredicate(pattern, target, mode))
        kwargs = {}
        if use_stereochemistry is not None:
            kwargs["do_chiral_sss"] = use_stereochemistry
        if similarity is not None:
            kwargs["tanimoto_threshold"] = similarity
        queries.append(ReactionComponentQuery(predicates, **kwargs))
    if not queries:
        raise ValueError("No query parameters were specified.")
    if limit:
        limit = min(limit, MAX_RESULTS)
    else:
        limit = MAX_RESULTS
    results = run_query(queries, limit)
    return QueryResult.from_results(results)


class ReactionQuery(BaseModel):
    """Reaction query."""

    reaction_ids: list[str]


@router.post("/api/fetch_reactions")
def fetch_reactions(inputs: ReactionQuery) -> list[QueryResult]:
    """Fetches a list of Reactions by ID."""
    command = ReactionIdQuery(inputs.reaction_ids)
    results = connect().run_query(command)
    return QueryResult.from_results(results)


class DatasetInfo(BaseModel):
    """Dataset info."""

    dataset_id: str
    name: str
    description: str
    size: int


@router.get("/api/fetch_datasets")
def fetch_datasets() -> list[DatasetInfo]:
    """Returns info about the current datasets."""
    with connect().connection as connection, connection.cursor() as cursor:
        cursor.execute(
            """
            SELECT dataset.dataset_id, name, description, size
            FROM dataset
            JOIN (
                SELECT dataset_id, COUNT(*) AS size
                FROM reaction
                GROUP BY dataset_id
            ) dataset_counts ON dataset.id = dataset_counts.dataset_id
            """
        )
        rows = []
        for row in cursor:
            rows.append(DatasetInfo(**row))
        return rows
