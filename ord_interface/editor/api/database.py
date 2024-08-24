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

"""Editor database."""
import os
from contextlib import contextmanager
from typing import Iterator

import psycopg
from ord_schema.orm.database import get_connection_string
from ord_schema.proto.dataset_pb2 import Dataset
from psycopg import Cursor
from psycopg.rows import dict_row

POSTGRES_DATABASE = "editor"


@contextmanager
def get_cursor() -> Iterator[Cursor]:
    """Returns a psycopg cursor."""
    dsn = os.getenv("ORD_EDITOR_POSTGRES")
    if dsn is None:
        dsn = get_connection_string(
            database=POSTGRES_DATABASE,
            username=os.environ["POSTGRES_USER"],
            password=os.environ["POSTGRES_PASSWORD"],
            host=os.environ["POSTGRES_HOST"],
        )
    with (  # pylint: disable=not-context-manager
        psycopg.connect(dsn, row_factory=dict_row) as connection,
        connection.cursor() as cursor,
    ):
        yield cursor


def prepare_database(cursor: Cursor) -> None:
    """Prepares the database."""
    cursor.execute(
        """
        CREATE TABLE users (
            user_id CHARACTER(32) PRIMARY KEY,
            user_name TEXT,
            created timestamp default current_timestamp NOT NULL
        )
        """
    )
    cursor.execute(
        """
        CREATE TABLE datasets (
            user_id CHARACTER(32) REFERENCES users,
            dataset_name TEXT NOT NULL,
            binpb BYTEA,
            created timestamp default current_timestamp NOT NULL,
            modified timestamp default current_timestamp NOT NULL,
            PRIMARY KEY (user_id, dataset_name)
        )
        """
    )
    # See https://aviyadav231.medium.com/98766e3b47a0.
    cursor.execute(
        """
        CREATE FUNCTION update_modified()
        RETURNS TRIGGER AS $$
        BEGIN
            NEW.modified = now();
            RETURN NEW;
        END;
        $$ language 'plpgsql'
        """
    )
    cursor.execute(
        """
        CREATE TRIGGER run_update_modified
            BEFORE UPDATE ON datasets FOR EACH ROW
        EXECUTE PROCEDURE update_modified();
        """
    )


def add_user(user_id: str, user_name: str, cursor: Cursor) -> None:
    """Adds a user to the database."""
    cursor.execute(
        """
        INSERT INTO users (user_id, user_name) VALUES (%s, %s)
            ON CONFLICT (user_id) DO UPDATE SET user_name = %s
        """,
        (user_id, user_name, user_name),
    )


def add_dataset(user_id: str, dataset: Dataset, cursor: Cursor) -> None:
    """Adds a dataset to the database."""
    binpb = dataset.SerializeToString()
    cursor.execute(
        """
        INSERT INTO datasets (user_id, dataset_name, binpb) VALUES (%s, %s, %s) 
            ON CONFLICT (user_id, dataset_name) DO UPDATE SET binpb = %s
        """,
        (user_id, dataset.name, binpb, binpb),
    )


def get_dataset(user_id: str, dataset_name: str, cursor: Cursor) -> Dataset | None:
    """Returns the requested dataset or None if the dataset does not exist."""
    cursor.execute("SELECT binpb FROM datasets WHERE user_id = %s AND dataset_name = %s", (user_id, dataset_name))
    row = cursor.fetchone()
    if row is None:
        return None
    return Dataset.FromString(row["binpb"])
