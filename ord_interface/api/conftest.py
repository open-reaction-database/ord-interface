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
from typing import Iterator

import psycopg2
import pytest
from psycopg2.extras import DictCursor
from sqlalchemy import create_engine
from sqlalchemy.orm import Session
from testing.postgresql import Postgresql

from ord_interface.api.testing import setup_test_postgres


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
def test_cursor(test_postgres) -> Iterator[DictCursor]:
    options = "-c search_path=public,ord"
    with psycopg2.connect(test_postgres.url(), cursor_factory=DictCursor, options=options) as connection:
        connection.set_session(readonly=True)
        with connection.cursor() as cursor:
            assert isinstance(cursor, DictCursor)  # Type hint.
            yield cursor
