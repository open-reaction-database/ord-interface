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

"""Open Reaction Database API."""

import gzip
import os
from base64 import b64decode, b64encode
from contextlib import asynccontextmanager, contextmanager
from io import BytesIO
from typing import Iterator
from uuid import uuid4

import psycopg
from fastapi import FastAPI, Request, Response, UploadFile
from google.protobuf import json_format, text_format  # pytype: disable=import-error
from ord_schema import resolvers
from ord_schema.message_helpers import create_message, mol_from_compound
from ord_schema.validations import ValidationOptions, validate_message
from ord_schema.proto.dataset_pb2 import Dataset
from ord_schema.proto.reaction_pb2 import Compound, Reaction
from ord_schema.orm.database import get_connection_string
from ord_schema.templating import generate_dataset, read_spreadsheet
from psycopg import Cursor
from psycopg.rows import dict_row
from testing.postgresql import Postgresql

from ord_interface.editor.api.database import add_dataset, get_dataset
from ord_interface.editor.api.testing import setup_test_postgres
from ord_interface.visualization.drawing import mol_to_svg
from ord_interface.visualization.generate_text import generate_html

POSTGRES_DATABASE = "editor"


@asynccontextmanager
async def lifespan(*args, **kwargs):
    """FastAPI lifespan setup; see https://fastapi.tiangolo.com/advanced/events/#lifespan."""
    del args, kwargs  # Unused.
    if os.getenv("ORD_EDITOR_TESTING", "FALSE") == "TRUE":
        with Postgresql() as postgres:
            setup_test_postgres(postgres.url())
            os.environ["ORD_EDITOR_POSTGRES"] = postgres.url()
            yield
    else:
        yield


app = FastAPI(lifespan=lifespan, root_path="/editor")


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
    with psycopg.connect(dsn, row_factory=dict_row) as connection, connection.cursor() as cursor:
        yield cursor


@app.get("/healthcheck")
async def health_check():
    return True


@app.get("/create_user", tags=["users"])
async def create_user(user_name: str):
    """Creates a new user and returns the associated user ID."""
    user_id = uuid4().hex
    with get_cursor() as cursor:
        cursor.execute("INSERT INTO users (user_id, user_name) VALUES (%s, %s)", (user_id, user_name))
    return user_id


@app.get("/delete_user", tags=["users"])
async def delete_user(user_id: str):
    """Deletes a user."""
    with get_cursor() as cursor:
        cursor.execute("DELETE FROM users WHERE user_id = %s", (user_id,))


@app.get("/list_datasets", tags=["datasets"])
async def list_datasets(user_id: str):
    """Returns a list of all datasets associated with the given user."""
    datasets = []
    with get_cursor() as cursor:
        cursor.execute("SELECT dataset_name FROM datasets WHERE user_id = %s", (user_id,))
        for row in cursor:
            datasets.append(row[0])
    return datasets


@app.get("/list_reactions", tags=["reactions"])
async def list_reactions(user_id: str, dataset_name: str):
    """Fetches a list of reactions in a dataset."""
    with get_cursor() as cursor:
        dataset = get_dataset(user_id, dataset_name, cursor)
    if dataset is None:
        return Response(status_code=404)
    return [reaction.reaction_id for reaction in dataset.reactions]


@app.get("/fetch_dataset", tags=["datasets"])
async def fetch_dataset(user_id: str, dataset_name: str):
    """Returns a base64-encoded dataset proto."""
    with get_cursor() as cursor:
        cursor.execute("SELECT binpb FROM datasets WHERE user_id = %s AND dataset_name = %s", (user_id, dataset_name))
        return b64encode(cursor.fetchone()["binpb"]).decode()


def send_message(message) -> str:
    """Converts a protocol buffer message to a base64-encoded string."""
    return b64encode(message.SerializeToString()).decode()


@app.get("/fetch_reaction", tags=["reactions"])
async def fetch_reaction(user_id: str, dataset_name: str, index: int):
    """Returns a base64-encoded reaction proto."""
    with get_cursor() as cursor:
        dataset = get_dataset(user_id, dataset_name, cursor)
    return send_message(dataset.reactions[index])


def download_message(message: Dataset | Reaction, prefix: str, kind: str):
    """Downloads a dataset or reaction as a gzipped file."""
    match kind:
        case "json":
            data = json_format.MessageToJson(message)
            filename = f"{prefix}.json"
        case "binpb":
            data = message.SerializeToString()
            filename = f"{prefix}.binpb"
        case "txtpb":
            data = text_format.MessageToBytes(message)
            filename = f"{prefix}.txtpb"
        case _:
            raise ValueError(kind)
    headers = {"Content-Disposition": f'attachment; filename="{filename}.gz"'}
    return Response(gzip.compress(data), headers=headers, media_type="application/gzip")


@app.get("/download_reaction", tags=["reactions"])
async def download_reaction(user_id: str, dataset_name: str, index: int, kind: str):
    """Downloads a reaction."""
    with get_cursor() as cursor:
        dataset = get_dataset(user_id, dataset_name, cursor)
    if dataset is None:
        return Response(status_code=404)
    return download_message(dataset.reactions[index], f"{dataset_name}-{index}", kind=kind)


@app.get("/download_dataset", tags=["datasets"])
async def download_dataset(user_id: str, dataset_name: str, kind: str):
    """Downloads a dataset."""
    # NOTE(skearnes): See https://protobuf.dev/reference/protobuf/textformat-spec/#text-format-files for comments on
    # preferred file extensions.
    with get_cursor() as cursor:
        dataset = get_dataset(user_id, dataset_name, cursor)
    if dataset is None:
        return Response(status_code=404)
    return download_message(dataset, dataset_name, kind=kind)


@app.post("/upload_dataset", tags=["datasets"])
async def upload_dataset(user_id: str, file: UploadFile):
    """Uploads a dataset."""
    data = await file.read()
    if file.filename.endswith(".gz"):
        data = gzip.decompress(data)
    if ".json" in file.filename:
        dataset = json_format.Parse(data.decode(), Dataset())
    elif ".binpb" in file.filename:
        dataset = Dataset.FromString(data)
    elif ".txtpb" in file.filename:
        dataset = text_format.Parse(data.decode(), Dataset())
    else:
        raise ValueError(file.filename)
    with get_cursor() as cursor:
        add_dataset(user_id, dataset, cursor)


@app.get("/create_dataset", tags=["datasets"])
async def create_dataset(user_id: str, dataset_name: str):
    """Creates a new dataset."""
    with get_cursor() as cursor:
        if get_dataset(user_id, dataset_name, cursor) is not None:
            return Response(status_code=409)
        cursor.execute("INSERT INTO datasets (user_id, dataset_name) VALUES (%s, %s)", (user_id, dataset_name))


@app.get("/delete_dataset", tags=["datasets"])
async def delete_dataset(user_id: str, dataset_name: str):
    """Deletes a dataset."""
    with get_cursor() as cursor:
        cursor.execute("DELETE FROM datasets WHERE user_id = %s AND dataset_name = %s", (user_id, dataset_name))


def adjust_error(error: str) -> str:
    """Strips the message name from errors to make them more readable."""
    fields = error.split(":")
    location = ".".join(fields[0].strip().split(".")[1:])
    message = ":".join(fields[1:])
    if location:
        return f"{location}: {message.strip()}"
    return message.strip()


@app.post("/validate/{message_type}", tags=["utilities"])
async def validate(message_type: str, request: Request):
    """Validates a protocol buffer message."""
    message = create_message(message_type)
    message.ParseFromString(await request.body())
    if message == type(message)():  # Skip empty messages.
        return {"errors": [], "warnings": []}
    options = ValidationOptions(require_provenance=True)
    output = validate_message(message, raise_on_error=False, options=options)
    errors = list(map(adjust_error, output.errors))
    warnings = list(map(adjust_error, output.warnings))
    return {"errors": errors, "warnings": warnings}


@app.post("/resolve_input", tags=["utilities"])
async def resolve_input(input_string: str):
    try:
        return send_message(resolvers.resolve_input(input_string))
    except (ValueError, KeyError) as error:
        return Response(str(error), status_code=400)


@app.post("/resolve_compound", tags=["utilities"])
async def resolve_compound(identifier_type: str, compound_string: str):
    try:
        smiles, resolver = resolvers.name_resolve(identifier_type, compound_string)
        return {"smiles": resolvers.canonicalize_smiles(smiles), "resolver": resolver}
    except ValueError as error:
        return Response(str(error), status_code=400)


@app.post("/canonicalize_smiles", tags=["utilities"])
async def canonicalize_smiles(smiles: str):
    try:
        return resolvers.canonicalize_smiles(smiles)
    except ValueError as error:
        return Response(str(error), status_code=400)


@app.post("/render_reaction", tags=["visualization"])
async def render_reaction(request: Request):
    """Returns a block of HTML that contains a visual summary of the reaction."""
    reaction = Reaction.FromString(await request.body())
    if not (reaction.inputs or reaction.outcomes):
        return ""
    try:
        return generate_html(reaction)
    except (ValueError, KeyError) as error:
        return Response(str(error), status_code=400)


@app.post("/render_compound", tags=["visualization"])
async def render_compound(request: Request):
    """Returns an HTML-tagged SVG for the given Compound."""
    compound = Compound.FromString(await request.body())
    try:
        return mol_to_svg(mol_from_compound(compound))
    except ValueError as error:
        return Response(str(error), status_code=400)


@app.post("/enumerate_dataset", tags=["utilities"])
async def enumerate_dataset(user_id: str, spreadsheet_name: str, spreadsheet_data: str, template_string: str):
    """Creates a new dataset based on a template reaction and a spreadsheet.

    A new dataset is created from the template and spreadsheet using ord_schema.templating.generate_dataset.

    Args:
        user_id: User ID.
        spreadsheet_name: The original filename of the uploaded spreadsheet.
        spreadsheet_data: A base64-encoded string containing the contents of the spreadsheet.
        template_string: A string containing a text-formatted Reaction proto, i.e., the contents of a txtpb file.
    """
    try:
        basename, suffix = os.path.splitext(spreadsheet_name)
        spreadsheet_data = BytesIO(b64decode(spreadsheet_data))
        dataframe = read_spreadsheet(spreadsheet_data, suffix=suffix)
        name = f"{basename}_dataset"
        dataset = generate_dataset(
            name=name,
            description="Enumerated by the ORD editor.",
            template_string=template_string,
            df=dataframe,
            validate=False,
        )
        with get_cursor() as cursor:
            add_dataset(user_id, dataset, cursor)
    except Exception as error:  # pylint: disable=broad-except
        return Response(str(error), status_code=400)
