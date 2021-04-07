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
"""Populates a PostgreSQL database containing the ORD.

For simplicity, and to aid in debugging, we write tables to individual CSV
files and load them into PostgreSQL with the COPY command.
"""

import dataclasses
import glob
import itertools
import os
import sys
from typing import Iterable, List, Mapping, Union

from absl import app
from absl import flags
from absl import logging
import psycopg2
from psycopg2 import extras
from psycopg2 import sql

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

import ord_interface

FLAGS = flags.FLAGS
flags.DEFINE_string('input', None, 'Input pattern (glob).')
flags.DEFINE_boolean('overwrite', False, 'If True, overwrite existing tables.')
flags.DEFINE_boolean('downsample', False,
                     'Whether to downsample datasets for testing.')
# Connection parameters.
flags.DEFINE_string('host', 'localhost', 'PostgreSQL server host.')
flags.DEFINE_string('dbname', ord_interface.POSTGRES_DB, 'Database name.')
flags.DEFINE_string('user', ord_interface.POSTGRES_USER, 'Username.')
flags.DEFINE_string('password', ord_interface.POSTGRES_PASSWORD, 'Password.')
flags.DEFINE_integer('port', ord_interface.POSTGRES_PORT, 'Port.')

# Maximum number of reactions to keep when downsampling datasets for testing.
_TEST_DATASET_SIZE = 100
# Datasets to use for testing.
_TEST_DATASETS = [
    'ord_dataset-d319c2a22ecf4ce59db1a18ae71d529c',
    'ord_dataset-cbcc4048add7468e850b6ec42549c70d',
    'ord_dataset-b440f8c90b6343189093770060fc4098',
    'ord_dataset-33320f511ffb4f89905448c7a5153111',
    'ord_dataset-7d8f5fd922d4497d91cb81489b052746',
    'ord_dataset-46ff9a32d9e04016b9380b1b1ef949c3',
]


@dataclasses.dataclass(frozen=True)
class InsertValues:
    """Container for INSERT VALUES for a single Reaction."""
    reactions: Mapping[str, Union[str, bytes]]
    inputs: Iterable[Mapping[str, str]]
    outputs: Iterable[Mapping[str, Union[str, float]]]


def process_reaction(reaction: reaction_pb2.Reaction,
                     dataset_id: str) -> InsertValues:
    """Adds a Reaction to the database.

    Args:
        reaction: Reaction proto.
        cursor: psycopg2 cursor.
        dataset_id: Dataset ID.

    Returns:
        InsertValues instance.
    """
    reactions_values = _reactions_table(reaction=reaction,
                                        dataset_id=dataset_id)
    inputs_values = _inputs_table(reaction=reaction)
    outputs_values = _outputs_table(reaction=reaction)
    return InsertValues(reactions_values, inputs_values, outputs_values)


def _reactions_table(reaction: reaction_pb2.Reaction,
                     dataset_id: str) -> Mapping[str, Union[str, bytes, None]]:
    """Adds a Reaction to the 'reactions' table.

    Args:
        reaction: Reaction proto.
        dataset_id: Dataset ID.

    Returns:
        Dict mapping string column names to values.
    """
    values = {
        'dataset_id': dataset_id,
        'reaction_id': reaction.reaction_id,
        'serialized': reaction.SerializeToString().hex()
    }
    try:
        reaction_smiles = message_helpers.get_reaction_smiles(
            reaction, generate_if_missing=True)
        # Control for REACTION_CXSMILES.
        values['reaction_smiles'] = reaction_smiles.split()[0]
    except ValueError:
        values['reaction_smiles'] = None
    if reaction.provenance.doi:
        values['doi'] = reaction.provenance.doi
    else:
        values['doi'] = None
    return values


def _inputs_table(reaction: reaction_pb2.Reaction) -> List[Mapping[str, str]]:
    """Adds rows to the 'inputs' table.

    Args:
        reaction: Reaction proto.

    Returns:
        List of dicts mapping string column names to values.
    """
    values = []
    for key in sorted(reaction.inputs):
        reaction_input = reaction.inputs[key]
        for compound in reaction_input.components:
            row = {'reaction_id': reaction.reaction_id}
            try:
                row['smiles'] = message_helpers.smiles_from_compound(compound)
            except ValueError:
                continue
            values.append(row)
    return values


def _outputs_table(
    reaction: reaction_pb2.Reaction
) -> List[Mapping[str, Union[str, float, None]]]:
    """Adds rows to the 'outputs' table.

    Args:
        reaction: Reaction proto.

    Returns:
        List of dicts mapping string column names to values.
    """
    values = []
    for outcome in reaction.outcomes:
        for product in outcome.products:
            row = {'reaction_id': reaction.reaction_id}
            try:
                row['smiles'] = message_helpers.smiles_from_compound(product)
            except ValueError:
                continue
            row['yield'] = message_helpers.get_product_yield(product)
            values.append(row)
    return values


def create_database(cursor: psycopg2.extensions.cursor, overwrite: bool):
    """Initializes the Postgres database."""
    if overwrite:
        logging.info('Removing existing tables')
        for table in ord_interface.TABLES:
            command = sql.SQL('DROP TABLE IF EXISTS {}')
            cursor.execute(command.format(sql.Identifier(table)))
    cursor.execute(sql.SQL('CREATE EXTENSION IF NOT EXISTS rdkit'))
    cursor.execute(sql.SQL('CREATE EXTENSION IF NOT EXISTS tsm_system_rows'))
    cursor.execute(
        sql.SQL('CREATE SCHEMA {}').format(
            sql.Identifier(ord_interface.RDKIT_SCHEMA)))
    for table, columns in ord_interface.TABLES.items():
        dtypes = []
        for column, dtype in columns.items():
            if table == 'reactions' and column == 'reaction_id':
                component = sql.SQL('{} {} PRIMARY KEY')
            else:
                component = sql.SQL('{} {}')
            # NOTE(kearnes): sql.Identifier(dtype) does not work for the
            # 'double precision' type.
            dtypes.append(
                component.format(sql.Identifier(column), sql.SQL(dtype)))
        command = sql.Composed([
            sql.SQL('CREATE TABLE {} (').format(sql.Identifier(table)),
            sql.Composed(dtypes).join(', '),
            sql.SQL(')')
        ])
        logging.info('Running:\n%s', command.as_string(cursor))
        cursor.execute(command)


def _rdkit_reaction_smiles(cursor: psycopg2.extensions.cursor, table: str):
    """Adds RDKit cartridge tables for reaction SMILES.

    Creates a new table rdk.<table> with the following columns:
        * r: reaction objects loaded from reaction SMILES

    An index is also created for each column.

    Args:
        cursor: psycopg2 cursor.
        table: Table name.
    """
    cursor.execute(
        sql.SQL("""
        SELECT reaction_id,
               r
        INTO {} FROM (
            SELECT reaction_id, 
                   reaction_from_smiles(reaction_smiles::cstring) AS r
            FROM {}) tmp
        WHERE r IS NOT NULL""").format(
            sql.Identifier(ord_interface.RDKIT_SCHEMA, table),
            sql.Identifier(table)))
    cursor.execute(
        sql.SQL('CREATE INDEX {} ON {} USING gist(r)').format(
            sql.Identifier(f'{table}_r'),
            sql.Identifier(ord_interface.RDKIT_SCHEMA, table)))


def _rdkit_smiles(cursor: psycopg2.extensions.cursor, table: str):
    """Adds RDKit cartridge tables for molecule SMILES.

    Creates a new table rdk.<table> with the following columns:
        * m: mol objects loaded from SMILES
        * mfp2: Morgan fingerprints

    An index is also created for each column.

    Args:
        cursor: psycopg2 cursor.
        table: Table name.
    """
    cursor.execute(
        sql.SQL("""
        SELECT reaction_id,
               m,
               morganbv_fp(m) AS mfp2 
        INTO {} FROM (
            SELECT reaction_id, 
                   mol_from_smiles(smiles::cstring) AS m
            FROM {}) tmp
        WHERE m IS NOT NULL""").format(
            sql.Identifier(ord_interface.RDKIT_SCHEMA, table),
            sql.Identifier(table)))
    cursor.execute(
        sql.SQL('CREATE INDEX {} ON {} USING gist(m)').format(
            sql.Identifier(f'{table}_m'),
            sql.Identifier(ord_interface.RDKIT_SCHEMA, table)))
    cursor.execute(
        sql.SQL('CREATE INDEX {} ON {} USING gist(mfp2)').format(
            sql.Identifier(f'{table}_mfp2'),
            sql.Identifier(ord_interface.RDKIT_SCHEMA, table)))


def process_dataset(filename: str, cursor: psycopg2.extensions.cursor,
                    downsample: bool):
    """Processes a single Dataset."""
    dataset_id = os.path.basename(filename).split('.')[0]
    if downsample and dataset_id not in _TEST_DATASETS:
        logging.info('TESTING: Dataset is not in _TEST_DATASETS')
        return
    dataset = message_helpers.load_message(filename, dataset_pb2.Dataset)
    if downsample and len(dataset.reactions) > _TEST_DATASET_SIZE:
        # Downsample ord-data Datasets for testing.
        logging.info('TESTING: Downsampling from %d->%d reactions',
                     len(dataset.reactions), _TEST_DATASET_SIZE)
        reactions = dataset.reactions[:_TEST_DATASET_SIZE]
    else:
        reactions = dataset.reactions
    values = []
    for reaction in reactions:
        values.append(process_reaction(reaction, dataset_id=dataset.dataset_id))
    # Update reactions table.
    extras.execute_values(cursor,
                          'INSERT INTO reactions VALUES %s',
                          [value.reactions for value in values],
                          template="""(%(reaction_id)s,
                                       %(reaction_smiles)s,
                                       %(doi)s,
                                       %(dataset_id)s,
                                       %(serialized)s)""")
    # Update inputs table.
    extras.execute_values(cursor,
                          'INSERT INTO inputs VALUES %s',
                          itertools.chain.from_iterable(
                              [value.inputs for value in values]),
                          template='(%(reaction_id)s, %(smiles)s)')
    # Update outputs table.
    extras.execute_values(cursor,
                          'INSERT INTO outputs VALUES %s',
                          itertools.chain.from_iterable(
                              [value.outputs for value in values]),
                          template='(%(reaction_id)s, %(smiles)s, %(yield)s)')


def main(argv):
    del argv  # Only used by app.run().
    filenames = glob.glob(FLAGS.input)
    logging.info('Found %d datasets', len(filenames))
    if not filenames:
        sys.exit(1)
    connection = psycopg2.connect(dbname=FLAGS.dbname,
                                  user=FLAGS.user,
                                  password=FLAGS.password,
                                  host=FLAGS.host,
                                  port=FLAGS.port)
    with connection:
        with connection.cursor() as cursor:
            create_database(overwrite=FLAGS.overwrite, cursor=cursor)
            for filename in filenames:
                logging.info(filename)
                process_dataset(filename=filename,
                                cursor=cursor,
                                downsample=FLAGS.downsample)
            for table, columns in ord_interface.TABLES.items():
                logging.info('Adding RDKit cartridge functionality')
                if 'reaction_smiles' in columns:
                    _rdkit_reaction_smiles(cursor, table)
                elif 'smiles' in columns:
                    _rdkit_smiles(cursor, table)
            connection.commit()


if __name__ == '__main__':
    flags.mark_flag_as_required('input')
    app.run(main)
