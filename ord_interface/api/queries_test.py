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
    run_queries,
)


def test_dataset_id_query(test_cursor):
    dataset_ids = ["ord_dataset-89b083710e2d441aa0040c361d63359f"]
    query = DatasetIdQuery(dataset_ids)
    results = run_queries(test_cursor, query, limit=10)
    assert len(results) == 10


def test_reaction_id_query(test_cursor):
    reaction_ids = ["ord-3f67aa5592fd434d97a577988d3fd241"]
    query = ReactionIdQuery(reaction_ids)
    results = run_queries(test_cursor, query, return_ids=False)
    assert [result.reaction_id for result in results] == reaction_ids


def test_reaction_smarts_query(test_cursor):
    query = ReactionSmartsQuery("[#6]>>[#7]")
    results = run_queries(test_cursor, query, limit=10, return_ids=True)
    assert len(results) == 10


def test_reaction_conversion_query(test_cursor):
    query = ReactionConversionQuery(min_conversion=50, max_conversion=90)
    results = run_queries(test_cursor, query, return_ids=True)
    assert len(results) == 7


def test_reaction_yield_query(test_cursor):
    query = ReactionYieldQuery(min_yield=50, max_yield=90)
    results = run_queries(test_cursor, query, limit=10, return_ids=True)
    assert len(results) == 10


def test_doi_query(test_cursor):
    dois = ["10.1126/science.1255525"]
    query = DoiQuery(dois)
    results = run_queries(test_cursor, query, limit=10)
    assert len(results) == 10
    for result in results:
        assert result.reaction.provenance.doi in dois


def test_exact_query(test_cursor):
    query = ReactionComponentQuery(
        "[Br]C1=CC=C(C(C)=O)C=C1", ReactionComponentQuery.Target.INPUT, ReactionComponentQuery.MatchMode.EXACT
    )
    results = run_queries(test_cursor, query, limit=5, return_ids=True)
    assert len(results) == 5


def test_substructure_query(test_cursor):
    query = ReactionComponentQuery(
        "C", ReactionComponentQuery.Target.OUTPUT, ReactionComponentQuery.MatchMode.SUBSTRUCTURE
    )
    results = run_queries(test_cursor, query, limit=10, return_ids=True)
    assert len(results) == 10


def test_chiral_substructure_query(test_cursor):
    query = ReactionComponentQuery(
        "OC1CC(O)C1", ReactionComponentQuery.Target.INPUT, ReactionComponentQuery.MatchMode.SUBSTRUCTURE
    )
    results = run_queries(test_cursor, query, limit=10, return_ids=True)
    assert len(results) == 10
    query = ReactionComponentQuery(
        "O[C@H]1C[C@H](O)C1", ReactionComponentQuery.Target.INPUT, ReactionComponentQuery.MatchMode.SUBSTRUCTURE
    )
    results = run_queries(test_cursor, query, limit=10, return_ids=True)
    assert len(results) == 10
    query = ReactionComponentQuery(
        "O[C@H]1C[C@H](O)C1",
        ReactionComponentQuery.Target.INPUT,
        ReactionComponentQuery.MatchMode.SUBSTRUCTURE,
        use_chirality=True,
    )
    results = run_queries(test_cursor, query, limit=10, return_ids=True)
    assert not results
    query = ReactionComponentQuery(
        "O[C@@H]1C[C@H](O)C1",
        ReactionComponentQuery.Target.INPUT,
        ReactionComponentQuery.MatchMode.SUBSTRUCTURE,
        use_chirality=True,
    )
    results = run_queries(test_cursor, query, limit=10, return_ids=True)
    assert len(results) == 10


def test_smarts_query(test_cursor):
    query = ReactionComponentQuery("[#6]", ReactionComponentQuery.Target.INPUT, ReactionComponentQuery.MatchMode.SMARTS)
    results = run_queries(test_cursor, query, limit=10, return_ids=True)
    assert len(results) == 10


def test_similarity_query(test_cursor):
    kwargs = {
        "pattern": "CC=O",
        "target": ReactionComponentQuery.Target.INPUT,
        "match_mode": ReactionComponentQuery.MatchMode.SIMILAR,
    }
    query = ReactionComponentQuery(**kwargs, similarity_threshold=0.5)
    results = run_queries(test_cursor, query, limit=10, return_ids=True)
    assert not results
    query = ReactionComponentQuery(**kwargs, similarity_threshold=0.05)
    results = run_queries(test_cursor, query, limit=10, return_ids=True)
    assert len(results) == 10


def test_bad_smiles(test_cursor):
    with pytest.raises(ValueError, match="Cannot parse pattern"):
        ReactionComponentQuery(
            "invalid_smiles", ReactionComponentQuery.Target.INPUT, ReactionComponentQuery.MatchMode.SUBSTRUCTURE
        )
