# Copyright 2020 Open Reaction Database Project Authors
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
"""Tests for ord_interface.build_database."""
import os

import docopt
import psycopg2
import pytest

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

import ord_interface
from ord_interface.client import build_database


def connect(dbname):
    return psycopg2.connect(
        dbname=dbname,
        user=ord_interface.client.POSTGRES_USER,
        password=ord_interface.client.POSTGRES_PASSWORD,
        host="localhost",
        port=ord_interface.client.POSTGRES_PORT,
    )


@pytest.fixture
def dataset_filename(tmp_path) -> str:
    # Create a test database.
    connection = connect(ord_interface.client.POSTGRES_DB)
    connection.set_session(autocommit=True)
    with connection.cursor() as cursor:
        cursor.execute("CREATE DATABASE test;")
    connection.close()
    # Create a test dataset.
    reaction = reaction_pb2.Reaction()
    reaction.reaction_id = "test"
    reaction.identifiers.add(value="reaction", type="REACTION_SMILES")
    input1 = reaction.inputs["input1"]
    input1.components.add().identifiers.add(value="input1", type="SMILES")
    input2 = reaction.inputs["input2"]
    input2.components.add().identifiers.add(value="input2a", type="SMILES")
    input2.components.add().identifiers.add(value="input2b", type="SMILES")
    outcome = reaction.outcomes.add()
    product = outcome.products.add()
    product.measurements.add(type="YIELD", percentage={"value": 2.5})
    product.identifiers.add(value="product", type="SMILES")
    reaction.provenance.doi = "10.0000/test.foo"
    dataset = dataset_pb2.Dataset(dataset_id="test_dataset", reactions=[reaction])
    dataset_filename = (tmp_path / "test.pb").as_posix()
    message_helpers.write_message(dataset, dataset_filename)
    yield dataset_filename
    # Remove the test database.
    connection = connect(ord_interface.client.POSTGRES_DB)
    connection.set_session(autocommit=True)
    with connection.cursor() as cursor:
        cursor.execute("DROP DATABASE test;")
    connection.close()


def test_main(dataset_filename):
    input_pattern = os.path.join(os.path.dirname(dataset_filename), "*.pb")
    argv = ["--input", input_pattern, "--dbname", "test"]
    build_database.main(docopt.docopt(build_database.__doc__, argv))
    # Sanity checks.
    connection = connect("test")
    with connection:
        with connection.cursor() as cursor:
            cursor.execute("SELECT * from reactions LIMIT 1;")
            row = cursor.fetchone()
            assert len(row) == 5
            cursor.execute("SELECT * from inputs LIMIT 1;")
            row = cursor.fetchone()
            assert len(row) == 2
            cursor.execute("SELECT * from outputs LIMIT 1;")
            row = cursor.fetchone()
            assert len(row) == 3
    connection.close()
