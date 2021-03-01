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

from absl.testing import absltest
from absl.testing import flagsaver
import pandas as pd
import psycopg2

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

import ord_interface
from ord_interface import build_database


class BuildDatabaseTest(absltest.TestCase):

    def setUp(self):
        super().setUp()
        # NOTE(kearnes): ord-postgres is the hostname in docker-compose.
        self.host = 'ord-postgres'
        # Create a test database.
        with self._connect(ord_interface.POSTGRES_DB) as connection:
            connection.set_isolation_level(
                psycopg2.extensions.ISOLATION_LEVEL_AUTOCOMMIT)
            with connection.cursor() as cursor:
                cursor.execute('CREATE DATABASE test;')
        # Create a test dataset.
        self.test_subdirectory = self.create_tempdir()
        reaction = reaction_pb2.Reaction()
        reaction.reaction_id = 'test'
        reaction.identifiers.add(value='reaction', type='REACTION_SMILES')
        input1 = reaction.inputs['input1']
        input1.components.add().identifiers.add(value='input1', type='SMILES')
        input2 = reaction.inputs['input2']
        input2.components.add().identifiers.add(value='input2a', type='SMILES')
        input2.components.add().identifiers.add(value='input2b', type='SMILES')
        outcome = reaction.outcomes.add()
        product = outcome.products.add()
        product.measurements.add(type='YIELD', percentage={'value': 2.5})
        product.identifiers.add(value='product', type='SMILES')
        reaction.provenance.doi = '10.0000/test.foo'
        self.dataset = dataset_pb2.Dataset(dataset_id='test_dataset',
                                           reactions=[reaction])
        message_helpers.write_message(
            self.dataset, os.path.join(self.test_subdirectory, 'test.pb'))

    def tearDown(self):
        # Remove the test database.
        with self._connect(ord_interface.POSTGRES_DB) as connection:
            connection.set_isolation_level(
                psycopg2.extensions.ISOLATION_LEVEL_AUTOCOMMIT)
            with connection.cursor() as cursor:
                cursor.execute('DROP DATABASE test;')

    def _connect(self, dbname):
        return psycopg2.connect(dbname=dbname,
                                user=ord_interface.POSTGRES_USER,
                                password=ord_interface.POSTGRES_PASSWORD,
                                host=self.host,
                                port=ord_interface.POSTGRES_PORT)

    def test_main(self):
        input_pattern = os.path.join(self.test_subdirectory, '*.pb')
        with flagsaver.flagsaver(input=input_pattern,
                                 dbname='test',
                                 host=self.host):
            build_database.main(())
        # Sanity checks.
        with self._connect('test') as connection:
            with connection.cursor() as cursor:
                cursor.execute('SELECT * from reactions LIMIT 1;')
                row = cursor.fetchone()
                self.assertLen(row, 5)
                cursor.execute('SELECT * from inputs LIMIT 1;')
                row = cursor.fetchone()
                self.assertLen(row, 2)
                cursor.execute('SELECT * from outputs LIMIT 1;')
                row = cursor.fetchone()
                self.assertLen(row, 3)


if __name__ == '__main__':
    absltest.main()
