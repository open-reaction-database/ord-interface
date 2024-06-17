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
import os
from typing import Iterator
from unittest.mock import patch

import pytest
from fastapi.testclient import TestClient

from ord_interface.api.main import app


@pytest.fixture(name="client", scope="session")
def client_fixture(test_postgres) -> Iterator[TestClient]:
    with TestClient(app) as client, patch.dict(os.environ, {"ORD_INTERFACE_POSTGRES": test_postgres.url()}):
        yield client


def test_get_datasets(client):
    response = client.get("/api/datasets")
    dataset_info = response.json()
    assert len(dataset_info) == 3
