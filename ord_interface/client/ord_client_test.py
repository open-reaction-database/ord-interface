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
"""Tests for ord_interface.ord_client."""
import pytest

from ord_interface.client import ord_client


@pytest.fixture
def client() -> ord_client.OrdClient:
    yield ord_client.OrdClient(target="http://localhost:5001")


@pytest.mark.parametrize("dataset_id,expected_num_reactions", (("ord_dataset-d319c2a22ecf4ce59db1a18ae71d529c", 264),))
def test_fetch_dataset(client, dataset_id, expected_num_reactions):
    dataset = client.fetch_dataset(dataset_id)
    assert len(dataset.reactions) == expected_num_reactions


@pytest.mark.parametrize(
    "dataset_ids,expected_num_reactions", ((["ord_dataset-d319c2a22ecf4ce59db1a18ae71d529c"], [264]),)
)
def test_fetch_datasets(client, dataset_ids, expected_num_reactions):
    datasets = client.fetch_datasets(dataset_ids)
    assert len(dataset_ids) == len(expected_num_reactions) == len(datasets)
    for dataset, expected in zip(datasets, expected_num_reactions):
        assert len(dataset.reactions) == expected


@pytest.mark.parametrize("reaction_id,created_by", (("ord-cf0d04017ede4c8aab8a15119c53e57b", "Connor W. Coley"),))
def test_fetch_reaction(client, reaction_id, created_by):
    reaction = client.fetch_reaction(reaction_id)
    assert reaction.provenance.record_created.person.name == created_by


@pytest.mark.parametrize(
    "reaction_ids,created_by",
    (
        (
            ["ord-cf0d04017ede4c8aab8a15119c53e57b", "ord-153de2a96e07484e93d525d81f966789"],
            ["Connor W. Coley", "Connor W. Coley"],
        ),
    ),
)
def test_fetch_reactions(client, reaction_ids, created_by):
    reactions = client.fetch_reactions(reaction_ids)
    assert len(reaction_ids) == len(created_by) == len(reactions)
    for reaction, expected in zip(reactions, created_by):
        assert reaction.provenance.record_created.person.name == expected


def test_query_dataset_ids(client):
    results = client.query(dataset_ids=["ord_dataset-89b083710e2d441aa0040c361d63359f"])
    assert len(results) == 24


def test_query_reaction_ids(client):
    results = client.query(
        reaction_ids=[
            "ord-cf0d04017ede4c8aab8a15119c53e57b",
            "ord-153de2a96e07484e93d525d81f966789",
        ]
    )
    assert len(results) == 2


def test_query_reaction_smarts(client):
    results = client.query(reaction_smarts="[Br]C1=CC=C(C(C)=O)C=C1>CN(C)C=O>")
    assert len(results) == 9


def test_query_dois(client):
    doi = "10.1126/science.1255525"
    results = client.query(dois=[doi])
    assert len(results) == 24
    for result in results:
        assert result.reaction.provenance.doi == doi


def test_query_single_component(client):
    component = ord_client.ComponentQuery("[Br]C1=CC=C(C(C)=O)C=C1", target="input", mode="exact")
    results = client.query(components=[component])
    assert len(results) == 10


def test_query_multiple_components(client):
    component1 = ord_client.ComponentQuery("[Br]C1=CC=C(C(C)=O)C=C1", target="input", mode="exact")
    component2 = ord_client.ComponentQuery("CC(C)(C)OC(=O)N1CCCCC1C(=O)O", target="input", mode="exact")
    component3 = ord_client.ComponentQuery("CC(=O)C1=CC=C(C2CCCCN2C(=O)OC(C)(C)C)C=C1", target="output", mode="exact")
    results = client.query(components=[component1, component2, component3])
    assert len(results) == 1


@pytest.mark.skip("TODO(skearnes): Add a test once we have more chiral stuff.")
def test_query_stereochemistry(client):
    pass


def test_query_similarity(client):
    component = ord_client.ComponentQuery("[Br]C1=CC=C(C(C)=O)C=C1", target="input", mode="similar")
    results = client.query(components=[component], similarity=0.6)
    assert len(results) == 10
