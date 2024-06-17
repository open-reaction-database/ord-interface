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

"""Tests for ord_interface.api.client."""
import gzip
import os
from typing import Iterator
from unittest.mock import patch

import pytest
from fastapi.testclient import TestClient
from ord_schema.proto import dataset_pb2

from ord_interface.api.main import app


@pytest.fixture(name="client", scope="session")
def client_fixture(test_postgres) -> Iterator[TestClient]:
    with TestClient(app) as client, patch.dict(os.environ, {"ORD_INTERFACE_POSTGRES": test_postgres.url()}):
        yield client


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
    expected = """
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2990    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  1  3  1  0
M  END
"""
    assert response.json() == expected


def test_get_search_results(client):
    response = client.post("/api/download_search_results", json={"reaction_ids": ["ord-3f67aa5592fd434d97a577988d3fd241"]})
    dataset = dataset_pb2.Dataset.FromString(gzip.decompress(response.read()))
    assert dataset.name == "ORD Search Results"
    assert len(dataset.reactions) == 1
    assert dataset.reactions[0].reaction_id == "ord-3f67aa5592fd434d97a577988d3fd241"
