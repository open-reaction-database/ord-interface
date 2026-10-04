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
from collections.abc import AsyncIterator, Iterator
from contextlib import ExitStack
from typing import Any
from unittest.mock import patch

import psycopg
import pytest
import pytest_asyncio
from fastapi.testclient import TestClient
from ord_schema.logging import get_logger
from psycopg import AsyncCursor
from psycopg.rows import dict_row
from testing.postgresql import Postgresql
from testing.redis import RedisServer

from ord_interface.api.main import app
from ord_interface.api.testing import setup_test_postgres

logger = get_logger(__name__)


@pytest.fixture(name="test_postgres", scope="session")
def test_postgres_fixture() -> Iterator[Postgresql]:
    with Postgresql() as postgres:
        setup_test_postgres(postgres.url())
        yield postgres


@pytest_asyncio.fixture
async def test_cursor(test_postgres) -> AsyncIterator[AsyncCursor[dict[str, Any]]]:
    async with await psycopg.AsyncConnection[dict[str, Any]].connect(
        test_postgres.url(), row_factory=dict_row, options="-c search_path=public,ord"
    ) as connection:
        await connection.set_read_only(True)
        async with connection.cursor() as cursor:
            yield cursor


@pytest.fixture(name="test_valkey", scope="session")
def test_valkey_fixture() -> Iterator[RedisServer]:
    """Runs a throwaway valkey-server and points ``get_valkey()`` at it.

    Every ``VALKEY_*`` variable is overridden, so a shell pointed at another server
    cannot leak into the tests.

    Yields:
        The running server.
    """
    with RedisServer(redis_server="valkey-server") as valkey_server:
        dsn = valkey_server.dsn()
        environment = {
            "VALKEY_HOST": dsn["host"],
            "VALKEY_PORT": str(dsn["port"]),
            "VALKEY_SSL": "0",
        }
        with patch.dict(os.environ, environment):
            yield valkey_server


@pytest.fixture(scope="session")
def test_client(test_postgres, test_valkey) -> Iterator[TestClient]:
    with TestClient(app) as client, ExitStack() as stack:
        # NOTE(skearnes): Set ORD_INTERFACE_POSTGRES to use that database instead of a testing.postgresql instance.
        # To force the use of testing.postgresl, set ORD_INTERFACE_TESTING=TRUE.
        if os.environ.get(
            "ORD_INTERFACE_TESTING", "FALSE"
        ) == "FALSE" and not os.environ.get("ORD_INTERFACE_POSTGRES"):
            stack.enter_context(
                patch.dict(os.environ, {"ORD_INTERFACE_POSTGRES": test_postgres.url()})
            )
        logger.debug(f"ORD_INTERFACE_POSTGRES={os.environ['ORD_INTERFACE_POSTGRES']}")
        yield client
