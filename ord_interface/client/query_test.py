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
"""Tests for ord_interface.query."""
import numpy as np
import pytest

import ord_interface
from ord_interface.client import query


@pytest.fixture(scope="module")
def connection() -> query.OrdPostgres:
    yield query.OrdPostgres(
        dbname=ord_interface.client.POSTGRES_DB,
        user=ord_interface.client.POSTGRES_USER,
        password=ord_interface.client.POSTGRES_PASSWORD,
        host="localhost",
        port=ord_interface.client.POSTGRES_PORT,
    )


def test_random_sample_query(connection):
    command = query.RandomSampleQuery(16)
    results = connection.run_query(command, return_ids=True)
    assert len(results) == 16


def test_dataset_id_query(connection):
    dataset_ids = ["ord_dataset-89b083710e2d441aa0040c361d63359f"]
    command = query.DatasetIdQuery(dataset_ids)
    results = connection.run_query(command, limit=10, return_ids=True)
    assert len(results) == 10


def test_reaction_id_query(connection):
    reaction_ids = ["ord-cf0d04017ede4c8aab8a15119c53e57b"]
    command = query.ReactionIdQuery(reaction_ids)
    results = connection.run_query(command, limit=10, return_ids=False)
    assert [result.reaction_id for result in results] == reaction_ids


def test_reaction_smarts_query(connection):
    pattern = "[#6]>>[#7]"
    command = query.ReactionSmartsQuery(pattern)
    results = connection.run_query(command, limit=10, return_ids=True)
    assert len(results) == 10


def test_reaction_conversion_query(connection):
    command = query.ReactionConversionQuery(min_conversion=50, max_conversion=90)
    results = connection.run_query(command, limit=10, return_ids=True)
    assert len(results) == 10


def test_reaction_yield_query(connection):
    command = query.ReactionYieldQuery(min_yield=50, max_yield=90)
    results = connection.run_query(command, limit=10, return_ids=True)
    assert len(results) == 10


def test_doi_query(connection):
    dois = ["10.1126/science.1255525"]
    command = query.DoiQuery(dois)
    results = connection.run_query(command, limit=10)
    assert len(results) == 10
    for result in results:
        assert result.reaction.provenance.doi in dois


def test_exact_query(connection):
    pattern = "[Br]C1=CC=C(C(C)=O)C=C1"
    predicates = [
        query.ReactionComponentPredicate(
            pattern,
            target=query.ReactionComponentPredicate.Target.INPUT,
            mode=query.ReactionComponentPredicate.MatchMode.EXACT,
        )
    ]
    command = query.ReactionComponentQuery(predicates)
    results = connection.run_query(command, limit=5, return_ids=True)
    assert len(results) == 5
    # Check that we remove redundant reaction IDs.
    reaction_ids = [result.reaction_id for result in results]
    assert len(reaction_ids) == len(np.unique(reaction_ids))


def test_substructure_query(connection):
    pattern = "C"
    predicates = [
        query.ReactionComponentPredicate(
            pattern,
            target=query.ReactionComponentPredicate.Target.INPUT,
            mode=query.ReactionComponentPredicate.MatchMode.SUBSTRUCTURE,
        )
    ]
    command = query.ReactionComponentQuery(predicates)
    results = connection.run_query(command, limit=10, return_ids=True)
    assert len(results) == 10
    # Check that we remove redundant reaction IDs.
    reaction_ids = [result.reaction_id for result in results]
    assert len(reaction_ids) == len(np.unique(reaction_ids))


def test_smarts_query(connection):
    pattern = "[#6]"
    predicates = [
        query.ReactionComponentPredicate(
            pattern,
            target=query.ReactionComponentPredicate.Target.INPUT,
            mode=query.ReactionComponentPredicate.MatchMode.SMARTS,
        )
    ]
    command = query.ReactionComponentQuery(predicates)
    results = connection.run_query(command, limit=10, return_ids=True)
    assert len(results) == 10
    # Check that we remove redundant reaction IDs.
    reaction_ids = [result.reaction_id for result in results]
    assert len(reaction_ids) == len(np.unique(reaction_ids))


def test_similarity_query(connection):
    pattern = "CC=O"
    predicates = [
        query.ReactionComponentPredicate(
            pattern,
            query.ReactionComponentPredicate.MatchMode.SMARTS,
            mode=query.ReactionComponentPredicate.MatchMode.SIMILAR,
        )
    ]
    command = query.ReactionComponentQuery(predicates, tanimoto_threshold=0.5)
    results = connection.run_query(command, limit=10, return_ids=True)
    assert not results
    command = query.ReactionComponentQuery(predicates, tanimoto_threshold=0.05)
    results = connection.run_query(command, limit=10, return_ids=True)
    assert len(results) == 10
    # Check that we remove redundant reaction IDs.
    reaction_ids = [result.reaction_id for result in results]
    assert len(reaction_ids) == len(np.unique(reaction_ids))


def test_bad_smiles(connection):
    pattern = "invalid_smiles"
    predicates = [
        query.ReactionComponentPredicate(
            pattern,
            query.ReactionComponentPredicate.MatchMode.SIMILAR,
            mode=query.ReactionComponentPredicate.MatchMode.SUBSTRUCTURE,
        )
    ]
    command = query.ReactionComponentQuery(predicates)
    with pytest.raises(query.QueryException, match="cannot parse pattern: invalid_smiles"):
        connection.run_query(command)
