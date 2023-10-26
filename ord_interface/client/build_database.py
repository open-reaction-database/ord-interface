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
"""Populates a PostgreSQL database containing the ORD (for local testing).

NOTE: For production, use ord_schema/orm/add_datasets.py instead of this script.

Usage:
    build_database.py --input=<str> [options]

Options:
    --input=<str>       Input pattern (glob)
    --host=<str>        PostgreSQL server host [default: localhost]
    --dbname=<str>      Database name
    --user=<str>        Username
    --password=<str>    Password
    --port=<int>        Port
"""
import glob
import logging

import docopt
import sqlalchemy
from sqlalchemy import orm

from ord_schema import message_helpers
from ord_schema.orm import database
from ord_schema.proto import dataset_pb2

import ord_interface.client

logger = logging.getLogger()


def main(kwargs):
    filenames = glob.glob(kwargs["--input"])
    logger.info("Found %d datasets", len(filenames))
    if not filenames:
        raise ValueError("--input did not match any files")
    connection_string = database.get_connection_string(
        database=kwargs["--dbname"] or ord_interface.client.POSTGRES_DB,
        username=kwargs["--user"] or ord_interface.client.POSTGRES_USER,
        password=kwargs["--password"] or ord_interface.client.POSTGRES_PASSWORD,
        host=kwargs["--host"],
        port=kwargs["--port"] or ord_interface.client.POSTGRES_PORT,
    )
    engine = sqlalchemy.create_engine(connection_string, future=True)
    with engine.begin() as connection:  # pytype: disable=attribute-error
        connection.execute(sqlalchemy.text("CREATE EXTENSION IF NOT EXISTS tsm_system_rows"))
    database.prepare_database(engine)
    with orm.Session(engine) as session:
        for filename in filenames:
            dataset = message_helpers.load_message(filename, dataset_pb2.Dataset)
            database.add_dataset(dataset, session)
            session.flush()
            database.update_rdkit_tables(dataset.dataset_id, session=session)
            session.flush()
            database.update_rdkit_ids(dataset.dataset_id, session=session)
            session.commit()


if __name__ == "__main__":
    main(docopt.docopt(__doc__))
