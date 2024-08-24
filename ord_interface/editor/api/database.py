"""Editor database."""

from ord_schema.proto.dataset_pb2 import Dataset
from psycopg import Cursor


def prepare_database(cursor: Cursor) -> None:
    """Prepares the database."""
    cursor.execute(
        """
        CREATE TABLE users (
            user_id CHARACTER(32) PRIMARY KEY,
            name TEXT,
            created timestamp default current_timestamp NOT NULL
        )
        """
    )
    cursor.execute(
        """
        CREATE TABLE logins (
            access_token TEXT PRIMARY KEY,
            user_id CHARACTER(32) REFERENCES users,
            timestamp timestamp default current_timestamp NOT NULL
        )
        """
    )
    cursor.execute(
        """
        CREATE TABLE datasets (
            user_id CHARACTER(32) REFERENCES users,
            name TEXT NOT NULL,
            serialized BYTEA NOT NULL,
            created timestamp default current_timestamp NOT NULL,
            PRIMARY KEY (user_id, name)
        )
        """
    )


def add_user(user_id: str, name: str, cursor: Cursor) -> None:
    """Adds a user to the database."""
    cursor.execute(
        """
        INSERT INTO users (user_id, name) VALUES (%s, %s)
            ON CONFLICT (user_id) DO UPDATE SET name = %s
        """,
        (user_id, name, name),
    )


def add_dataset(user_id: str, dataset: Dataset, cursor: Cursor) -> None:
    """Adds a dataset to the database."""
    serialized = dataset.SerializeToString()
    cursor.execute(
        """
        INSERT INTO datasets (user_id, name, serialized) VALUES (%s, %s, %s) 
            ON CONFLICT (user_id, name) DO UPDATE SET serialized = %s
        """,
        (user_id, dataset.name, serialized, serialized),
    )
