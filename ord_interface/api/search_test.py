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
import os
from typing import Iterator
from unittest.mock import patch

import pytest
from fastapi.testclient import TestClient
from ord_schema.proto import dataset_pb2
from rdkit import Chem

from ord_interface.api.main import app


@pytest.fixture(name="client", scope="session")
def client_fixture(test_postgres) -> Iterator[TestClient]:
    with TestClient(app) as client, patch.dict(os.environ, {"ORD_INTERFACE_POSTGRES": test_postgres.url()}):
        yield client


@pytest.mark.parametrize(
    "params,num_expected",
    [
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
    ],
)
def test_query(client, params, num_expected):
    response = client.get("/api/query", params=params)
    assert len(response.json()) == num_expected


def test_get_reactions(client):
    response = client.post("/api/reactions", json={"reaction_ids": ["ord-3f67aa5592fd434d97a577988d3fd241"]})
    reactions = response.json()
    assert len(reactions) == 1
    assert reactions[0]["dataset_id"] == "ord_dataset-89b083710e2d441aa0040c361d63359f"


def test_get_datasets(client):
    response = client.get("/api/datasets")
    dataset_info = response.json()
    assert len(dataset_info) == 3


def test_get_molfile(client):
    response = client.get("/api/molfile", params={"smiles": "C(=O)N"})
    assert Chem.MolToSmiles(Chem.MolFromMolBlock(response.json())) == "NC=O"


def test_get_search_results(client):
    response = client.post(
        "/api/download_search_results", json={"reaction_ids": ["ord-3f67aa5592fd434d97a577988d3fd241"]}
    )
    dataset = dataset_pb2.Dataset.FromString(gzip.decompress(response.read()))
    assert dataset.name == "ORD Search Results"
    assert len(dataset.reactions) == 1
    assert dataset.reactions[0].reaction_id == "ord-3f67aa5592fd434d97a577988d3fd241"