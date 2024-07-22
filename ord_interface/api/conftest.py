# Copyright 2022 Open Reaction Database Project Authors
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

"""Pytest fixtures."""
import os
from contextlib import ExitStack
from typing import Iterator
from unittest.mock import patch

import psycopg
import pytest
from fastapi.testclient import TestClient
from ord_schema.logging import get_logger
from psycopg import Cursor
from psycopg.rows import dict_row
from sqlalchemy import create_engine
from sqlalchemy.orm import Session
from testing.postgresql import Postgresql

from ord_interface.api.main import app
from ord_interface.api.testing import setup_test_postgres

logger = get_logger(__name__)


@pytest.fixture(name="test_postgres", scope="session")
def test_postgres_fixture() -> Iterator[Postgresql]:
    with Postgresql() as postgres:
        setup_test_postgres(postgres.url())
        yield postgres


@pytest.fixture
def test_session(test_postgres) -> Iterator[Session]:
    engine = create_engine(test_postgres.url(), future=True)
    with Session(engine) as session:
        yield session


@pytest.fixture
def test_cursor(test_postgres) -> Iterator[Cursor]:
    options = "-c search_path=public,ord"
    with psycopg.connect(test_postgres.url(), row_factory=dict_row, options=options) as connection:
        connection.set_session(readonly=True)
        with connection.cursor() as cursor:
            yield cursor


@pytest.fixture(scope="session")
def test_client(test_postgres) -> Iterator[TestClient]:
    with TestClient(app) as client, ExitStack() as stack:
        # NOTE(skearnes): Set ORD_INTERFACE_POSTGRES to use that database instead of a testing.postgresql instance.
        # To force the use of testing.postgresl, set ORD_INTERFACE_TESTING=TRUE.
        if os.environ.get("ORD_INTERFACE_TESTING", "FALSE") == "FALSE" and not os.environ.get("ORD_INTERFACE_POSTGRES"):
            stack.enter_context(patch.dict(os.environ, {"ORD_INTERFACE_POSTGRES": test_postgres.url()}))
        logger.info(f"ORD_INTERFACE_POSTGRES={os.environ['ORD_INTERFACE_POSTGRES']}")
        yield client
