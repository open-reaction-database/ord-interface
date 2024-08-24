"""Editor database."""

from ord_schema.proto.dataset_pb2 import Dataset
from psycopg import Cursor


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
