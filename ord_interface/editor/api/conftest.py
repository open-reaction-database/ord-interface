# Copyright 2024 Open Reaction Database Project Authors
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
from testing.postgresql import Postgresql
from testing.redis import RedisServer

from ord_interface.editor.api.main import app
from ord_interface.editor.api.testing import setup_test_postgres

logger = get_logger(__name__)


@pytest.fixture(name="test_postgres", scope="session")
def test_postgres_fixture() -> Iterator[Postgresql]:
    with Postgresql() as postgres:
        setup_test_postgres(postgres.url())
        yield postgres


@pytest.fixture
def test_cursor(test_postgres) -> Iterator[Cursor]:
    with (  # pylint: disable=not-context-manager
        psycopg.connect(test_postgres.url(), row_factory=dict_row) as connection,
        connection.cursor() as cursor,
    ):
        yield cursor


@pytest.fixture(scope="session")
def test_client(test_postgres) -> Iterator[TestClient]:
    with (
        TestClient(app) as client,
        RedisServer() as redis_server,
        patch.dict(os.environ, {"REDIS_PORT": str(redis_server.dsn()["port"])}),
        ExitStack() as stack,
    ):
        # NOTE(skearnes): Set ORD_EDITOR_POSTGRES to use that database instead of a testing.postgresql instance.
        # To force the use of testing.postgresl, set ORD_EDITOR_TESTING=TRUE.
        if os.environ.get("ORD_EDITOR_TESTING", "FALSE") == "FALSE" and not os.environ.get("ORD_EDITOR_POSTGRES"):
            stack.enter_context(patch.dict(os.environ, {"ORD_EDITOR_POSTGRES": test_postgres.url()}))
        logger.debug(f"ORD_EDITOR_POSTGRES={os.environ['ORD_EDITOR_POSTGRES']}")
        yield client
