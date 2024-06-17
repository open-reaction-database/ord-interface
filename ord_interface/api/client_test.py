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
