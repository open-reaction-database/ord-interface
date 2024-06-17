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
