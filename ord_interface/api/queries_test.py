# Copyright 2020 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Tests for ord_interface.api.queries."""
import pytest

from ord_interface.api.queries import (
    DatasetIdQuery,
    DoiQuery,
    ReactionComponentQuery,
    ReactionConversionQuery,
    ReactionIdQuery,
    ReactionSmartsQuery,
    ReactionYieldQuery,
    fetch_reactions,
    run_queries,
)


@pytest.mark.asyncio
async def test_fetch_reactions(test_cursor):
    results = await fetch_reactions(test_cursor, ["ord-1e8382606d99485b9859da6a92f80a72"])
    assert len(results) == 1
    result = results[0]
    assert result.dataset_id == "ord_dataset-b440f8c90b6343189093770060fc4098"
    assert result.reaction_id == "ord-1e8382606d99485b9859da6a92f80a72"
    assert result.reaction.reaction_id == result.reaction_id


@pytest.mark.asyncio
async def test_dataset_id_query(test_cursor):
    dataset_ids = ["ord_dataset-89b083710e2d441aa0040c361d63359f"]
    query = DatasetIdQuery(dataset_ids)
    results = await run_queries(test_cursor, query, limit=10)
    assert len(results) == 10


@pytest.mark.asyncio
async def test_reaction_id_query(test_cursor):
    reaction_ids = ["ord-3f67aa5592fd434d97a577988d3fd241"]
    query = ReactionIdQuery(reaction_ids)
    results = await run_queries(test_cursor, query)
    assert results == reaction_ids


@pytest.mark.asyncio
async def test_reaction_smarts_query(test_cursor):
    query = ReactionSmartsQuery("[#6]>>[#7]")
    results = await run_queries(test_cursor, query, limit=10)
    assert len(results) == 10


@pytest.mark.asyncio
async def test_reaction_conversion_query(test_cursor):
    query = ReactionConversionQuery(min_conversion=50, max_conversion=90)
    results = await run_queries(test_cursor, query)
    assert len(results) == 7


@pytest.mark.asyncio
async def test_reaction_yield_query(test_cursor):
    query = ReactionYieldQuery(min_yield=50, max_yield=90)
    results = await run_queries(test_cursor, query, limit=10)
    assert len(results) == 10


@pytest.mark.asyncio
async def test_doi_query(test_cursor):
    dois = ["10.1126/science.1255525"]
    query = DoiQuery(dois)
    results = await fetch_reactions(test_cursor, await run_queries(test_cursor, query, limit=10))
    assert len(results) == 10
    for result in results:
        assert result.reaction.provenance.doi in dois


@pytest.mark.asyncio
async def test_exact_query(test_cursor):
    query = ReactionComponentQuery(
        "[Br]C1=CC=C(C(C)=O)C=C1", ReactionComponentQuery.Target.INPUT, ReactionComponentQuery.MatchMode.EXACT
    )
    results = await run_queries(test_cursor, query, limit=5)
    assert len(results) == 5


@pytest.mark.asyncio
async def test_substructure_query(test_cursor):
    query = ReactionComponentQuery(
        "C", ReactionComponentQuery.Target.OUTPUT, ReactionComponentQuery.MatchMode.SUBSTRUCTURE
    )
    results = await run_queries(test_cursor, query, limit=10)
    assert len(results) == 10


@pytest.mark.asyncio
async def test_chiral_substructure_query(test_cursor):
    query = ReactionComponentQuery(
        "OC1CC(O)C1", ReactionComponentQuery.Target.INPUT, ReactionComponentQuery.MatchMode.SUBSTRUCTURE
    )
    results = await run_queries(test_cursor, query, limit=10)
    assert len(results) == 10
    query = ReactionComponentQuery(
        "O[C@H]1C[C@H](O)C1", ReactionComponentQuery.Target.INPUT, ReactionComponentQuery.MatchMode.SUBSTRUCTURE
    )
    results = await run_queries(test_cursor, query, limit=10)
    assert len(results) == 10
    query = ReactionComponentQuery(
        "O[C@H]1C[C@H](O)C1",
        ReactionComponentQuery.Target.INPUT,
        ReactionComponentQuery.MatchMode.SUBSTRUCTURE,
        use_chirality=True,
    )
    results = await run_queries(test_cursor, query, limit=10)
    assert not results
    query = ReactionComponentQuery(
        "O[C@@H]1C[C@H](O)C1",
        ReactionComponentQuery.Target.INPUT,
        ReactionComponentQuery.MatchMode.SUBSTRUCTURE,
        use_chirality=True,
    )
    results = await run_queries(test_cursor, query, limit=10)
    assert len(results) == 10


@pytest.mark.asyncio
async def test_smarts_query(test_cursor):
    query = ReactionComponentQuery("[#6]", ReactionComponentQuery.Target.INPUT, ReactionComponentQuery.MatchMode.SMARTS)
    results = await run_queries(test_cursor, query, limit=10)
    assert len(results) == 10


@pytest.mark.asyncio
async def test_similarity_query(test_cursor):
    kwargs = {
        "pattern": "CC=O",
        "target": ReactionComponentQuery.Target.INPUT,
        "match_mode": ReactionComponentQuery.MatchMode.SIMILAR,
    }
    query = ReactionComponentQuery(**kwargs, similarity_threshold=0.5)
    results = await run_queries(test_cursor, query, limit=10)
    assert not results
    query = ReactionComponentQuery(**kwargs, similarity_threshold=0.05)
    results = await run_queries(test_cursor, query, limit=10)
    assert len(results) == 10


@pytest.mark.asyncio
async def test_bad_smiles(test_cursor):
    with pytest.raises(ValueError, match="Cannot parse pattern"):
        ReactionComponentQuery(
            "invalid_smiles", ReactionComponentQuery.Target.INPUT, ReactionComponentQuery.MatchMode.SUBSTRUCTURE
        )
