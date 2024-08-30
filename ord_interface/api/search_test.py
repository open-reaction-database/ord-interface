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

"""Tests for ord_interface.api.search."""
import gzip

import pytest
from ord_schema.proto import dataset_pb2
from rdkit import Chem
from tenacity import retry, stop_after_attempt, wait_fixed

from ord_interface.api.queries import QueryResult

QUERY_PARAMS = [
    # Single factor queries.
    ({"dataset_id": ["ord_dataset-89b083710e2d441aa0040c361d63359f"]}, 24),
    ({"reaction_id": ["ord-3f67aa5592fd434d97a577988d3fd241"]}, 1),
    ({"reaction_smarts": "[#6]>>[#7]"}, 83),
    ({"min_conversion": 50, "max_conversion": 90}, 7),
    ({"min_yield": 50, "max_yield": 90}, 51),
    ({"doi": ["10.1126/science.1255525"]}, 24),
    ({"component": ["[Br]C1=CC=C(C(C)=O)C=C1;input;exact"]}, 10),
    ({"component": ["C;input;substructure"]}, 144),
    ({"component": ["O[C@@H]1C[C@H](O)C1;input;substructure"], "use_stereochemistry": True}, 20),
    ({"component": ["[#6];input;smarts"]}, 144),
    ({"component": ["CC=O;input;similar"], "similarity": 0.5}, 0),
    ({"component": ["CC=O;input;similar"], "similarity": 0.05}, 120),
    # Multi-factor queries.
    (
        {
            "min_yield": 50,
            "max_yield": 90,
            "component": ["[Br]C1=CC=C(C(C)=O)C=C1;input;exact", "CC(C)(C)OC(=O)NC;input;substructure"],
        },
        7,
    ),
]


@pytest.mark.parametrize("params,num_expected", QUERY_PARAMS)
def test_query(test_client, params, num_expected):
    response = test_client.get("/api/query", params=params)
    response.raise_for_status()
    assert len(response.json()) == num_expected


def test_get_reaction(test_client):
    response = test_client.get("/api/reaction", params={"reaction_id": "ord-3f67aa5592fd434d97a577988d3fd241"})
    response.raise_for_status()
    result = QueryResult(**response.json())
    assert result.dataset_id == "ord_dataset-89b083710e2d441aa0040c361d63359f"
    assert result.reaction.reaction_id == result.reaction_id


def test_get_reactions(test_client):
    response = test_client.post("/api/reactions", json={"reaction_ids": ["ord-3f67aa5592fd434d97a577988d3fd241"]})
    response.raise_for_status()
    reactions = response.json()
    assert len(reactions) == 1
    assert reactions[0]["dataset_id"] == "ord_dataset-89b083710e2d441aa0040c361d63359f"


def test_get_datasets(test_client):
    response = test_client.get("/api/datasets")
    response.raise_for_status()
    dataset_info = response.json()
    assert len(dataset_info) == 3


def test_get_molfile(test_client):
    response = test_client.get("/api/molfile", params={"smiles": "C(=O)N"})
    response.raise_for_status()
    assert Chem.MolToSmiles(Chem.MolFromMolBlock(response.json())) == "NC=O"


def test_get_search_results(test_client):
    response = test_client.post(
        "/api/download_search_results", json={"reaction_ids": ["ord-3f67aa5592fd434d97a577988d3fd241"]}
    )
    response.raise_for_status()
    dataset = dataset_pb2.Dataset.FromString(gzip.decompress(response.read()))
    assert dataset.name == "ORD Search Results"
    assert len(dataset.reactions) == 1
    assert dataset.reactions[0].reaction_id == "ord-3f67aa5592fd434d97a577988d3fd241"


@retry(stop=stop_after_attempt(5), wait=wait_fixed(1))
def wait_for_task(client, task_id) -> list[QueryResult]:
    response = client.get("/api/fetch_query_result", params={"task_id": task_id})
    response.raise_for_status()
    return response.json()


@pytest.mark.parametrize("params,num_expected", QUERY_PARAMS)
def test_query_async(test_client, params, num_expected):
    response = test_client.get("/api/submit_query", params=params)
    response.raise_for_status()
    assert len(wait_for_task(test_client, response.json())) == num_expected
