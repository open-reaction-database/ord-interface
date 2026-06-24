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

from typing import Any

import pytest

from ord_interface.api.queries import (
    DatasetIdQuery,
    DoiQuery,
    ReactionComponentQuery,
    ReactionConversionQuery,
    ReactionIdQuery,
    ReactionSmartsQuery,
    ReactionYieldQuery,
    fetch_dataset_most_used_smiles_for_inputs,
    fetch_dataset_most_used_smiles_for_products,
    fetch_reactions,
    run_queries,
)


@pytest.mark.asyncio
async def test_fetch_reactions(test_cursor):
    results = await fetch_reactions(
        test_cursor, ["ord-1e8382606d99485b9859da6a92f80a72"]
    )
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
    results = await fetch_reactions(
        test_cursor, await run_queries(test_cursor, query, limit=10)
    )
    assert len(results) == 10
    for result in results:
        assert result.reaction.provenance.doi in dois


@pytest.mark.asyncio
async def test_exact_query(test_cursor):
    query = ReactionComponentQuery(
        "[Br]C1=CC=C(C(C)=O)C=C1",
        ReactionComponentQuery.Target.INPUT,
        ReactionComponentQuery.MatchMode.EXACT,
    )
    results = await run_queries(test_cursor, query, limit=5)
    assert len(results) == 5


@pytest.mark.asyncio
async def test_substructure_query(test_cursor):
    query = ReactionComponentQuery(
        "C",
        ReactionComponentQuery.Target.OUTPUT,
        ReactionComponentQuery.MatchMode.SUBSTRUCTURE,
    )
    results = await run_queries(test_cursor, query, limit=10)
    assert len(results) == 10


@pytest.mark.asyncio
async def test_chiral_substructure_query(test_cursor):
    query = ReactionComponentQuery(
        "OC1CC(O)C1",
        ReactionComponentQuery.Target.INPUT,
        ReactionComponentQuery.MatchMode.SUBSTRUCTURE,
    )
    results = await run_queries(test_cursor, query, limit=10)
    assert len(results) == 10
    query = ReactionComponentQuery(
        "O[C@H]1C[C@H](O)C1",
        ReactionComponentQuery.Target.INPUT,
        ReactionComponentQuery.MatchMode.SUBSTRUCTURE,
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
    query = ReactionComponentQuery(
        "[#6]",
        ReactionComponentQuery.Target.INPUT,
        ReactionComponentQuery.MatchMode.SMARTS,
    )
    results = await run_queries(test_cursor, query, limit=10)
    assert len(results) == 10


@pytest.mark.asyncio
async def test_similarity_query(test_cursor):
    kwargs: dict[str, Any] = {
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


async def _max_input_similarities(test_cursor, reaction_ids, pattern):
    """Returns each reaction's greatest input-component Tanimoto similarity to ``pattern``.

    Independently recomputes the ranking score for INPUT-target similarity queries
    (it hardcodes the input-component join); an OUTPUT-target check would need the
    product-component join instead.
    """
    await test_cursor.execute(
        """
        SELECT reaction.reaction_id,
               MAX(tanimoto_sml(rdkit.mols.morgan_bfp, morganbv_fp(%s))) AS similarity
        FROM ord.reaction
        JOIN ord.reaction_input ON reaction_input.reaction_id = reaction.id
        JOIN ord.compound ON compound.reaction_input_id = reaction_input.id
        JOIN rdkit.mols ON rdkit.mols.id = compound.rdkit_mol_id
        WHERE reaction.reaction_id = ANY (%s)
        GROUP BY reaction.reaction_id
        """,
        [pattern, reaction_ids],
    )
    return {row["reaction_id"]: row["similarity"] async for row in test_cursor}


@pytest.mark.asyncio
async def test_similarity_ranking(test_cursor):
    kwargs: dict[str, Any] = {
        "pattern": "CC=O",
        "target": ReactionComponentQuery.Target.INPUT,
        "match_mode": ReactionComponentQuery.MatchMode.SIMILAR,
    }
    query = ReactionComponentQuery(**kwargs, similarity_threshold=0.05)
    ranked = await run_queries(test_cursor, query)
    assert len(ranked) > 10  # Enough matches to make the top-N check meaningful.
    # Results are ordered by descending best-component similarity.
    scores = await _max_input_similarities(test_cursor, ranked, "CC=O")
    ordered = [scores[reaction_id] for reaction_id in ranked]
    assert ordered == sorted(ordered, reverse=True)
    # Limiting returns the most similar matches, not an arbitrary subset.
    top = await run_queries(test_cursor, query, limit=10)
    assert top == ranked[:10]


@pytest.mark.asyncio
async def test_bad_smiles(test_cursor):
    with pytest.raises(ValueError, match="Cannot parse pattern"):
        ReactionComponentQuery(
            "invalid_smiles",
            ReactionComponentQuery.Target.INPUT,
            ReactionComponentQuery.MatchMode.SUBSTRUCTURE,
        )


@pytest.mark.asyncio
async def test_fetch_dataset_most_used_smiles_for_inputs(test_cursor):
    dataset_id = "ord_dataset-89b083710e2d441aa0040c361d63359f"
    results = await fetch_dataset_most_used_smiles_for_inputs(
        test_cursor, dataset_id, limit=10
    )
    assert len(results) == 10


@pytest.mark.asyncio
async def test_fetch_dataset_most_used_smiles_for_products(test_cursor):
    dataset_id = "ord_dataset-89b083710e2d441aa0040c361d63359f"
    results = await fetch_dataset_most_used_smiles_for_products(
        test_cursor, dataset_id, limit=10
    )
    assert len(results) == 10
