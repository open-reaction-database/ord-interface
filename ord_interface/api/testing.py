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

"""Testing utilities."""
import os
from glob import glob

from ord_schema.message_helpers import load_message
from ord_schema.orm.database import add_dataset, prepare_database, update_rdkit_ids, update_rdkit_tables
from ord_schema.proto import dataset_pb2
from sqlalchemy import create_engine
from sqlalchemy.orm import Session


def setup_test_postgres(url: str) -> None:
    datasets = [
        load_message(filename, dataset_pb2.Dataset)
        for filename in glob(os.path.join(os.path.dirname(__file__), "testdata", "*.pb.gz"))
    ]
    engine = create_engine(url, future=True)
    rdkit_cartridge = prepare_database(engine)
    with Session(engine) as session:
        for dataset in datasets:
            add_dataset(dataset, session)
            if rdkit_cartridge:
                session.flush()
                update_rdkit_tables(dataset.dataset_id, session)
                session.flush()
                update_rdkit_ids(dataset.dataset_id, session)
            session.commit()
