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
from contextlib import contextmanager
from typing import Annotated, Iterator

import psycopg2
from fastapi import APIRouter, Query, Response
from ord_schema.proto import dataset_pb2
from psycopg2.extras import DictCursor
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
def get_cursor() -> Iterator[DictCursor]:
    """Returns a psycopg2 cursor."""
    kwargs = {
        "dsn": os.getenv("ORD_INTERFACE_POSTGRES", "postgresql://postgres@localhost:5432/ord"),
        "cursor_factory": DictCursor,
        "options": "-c search_path=public,ord",
    }
    with psycopg2.connect(**kwargs) as connection:
        connection.set_session(readonly=True)
        with connection.cursor() as cursor:
            yield cursor


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
        kwargs = {}
        if use_stereochemistry is not None:
            kwargs["use_chirality"] = use_stereochemistry
        if similarity is not None:
            kwargs["similarity_threshold"] = similarity
        for spec in component:
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
    if limit:
        limit = min(limit, MAX_RESULTS)
    else:
        limit = MAX_RESULTS
    with get_cursor() as cursor:
        results = run_queries(cursor, queries, limit=limit)
    return QueryResult.from_results(results)


class ReactionIdList(BaseModel):
    """Reaction ID input."""

    reaction_ids: list[str]


@router.post("/reactions")
def get_reactions(inputs: ReactionIdList) -> list[QueryResult]:
    """Fetches a list of Reactions by ID."""
    with get_cursor() as cursor:
        results = run_queries(cursor, ReactionIdQuery(inputs.reaction_ids))
    return QueryResult.from_results(results)


class DatasetInfo(BaseModel):
    """Dataset info."""

    dataset_id: str
    name: str
    description: str | None
    num_reactions: int


@router.get("/datasets")
def get_datasets() -> list[DatasetInfo]:
    """Returns info about the current datasets."""
    with get_cursor() as cursor:
        cursor.execute("SELECT dataset_id, name, description, num_reactions FROM dataset")
        return [DatasetInfo(**row) for row in cursor]


@router.get("/molfile")
def get_molfile(smiles: str) -> str:
    """Returns a molblock for the given SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(smiles)
    return Chem.MolToMolBlock(mol)


@router.post("/search_results")
def get_search_results(inputs: ReactionIdList):
    """Downloads search results as a Dataset proto."""
    results = get_reactions(inputs)
    dataset = dataset_pb2.Dataset(name="ORD Search Results", reactions=[result.reaction for result in results])
    return Response(gzip.compress(dataset.SerializeToString()), media_type="application/gzip")
