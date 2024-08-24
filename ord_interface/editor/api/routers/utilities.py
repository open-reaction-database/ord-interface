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

"""Utility API endpoints."""

import os
from base64 import b64decode
from io import BytesIO

from fastapi import APIRouter, Request, Response
from ord_schema import resolvers
from ord_schema.message_helpers import create_message
from ord_schema.templating import generate_dataset, read_spreadsheet
from ord_schema.validations import ValidationOptions, validate_message

from ord_interface.editor.api import send_message
from ord_interface.editor.api.database import add_dataset, get_cursor

router = APIRouter(tags=["utilities"])


def adjust_error(error: str) -> str:
    """Strips the message name from errors to make them more readable."""
    fields = error.split(":")
    location = ".".join(fields[0].strip().split(".")[1:])
    message = ":".join(fields[1:])
    if location:
        return f"{location}: {message.strip()}"
    return message.strip()


@router.post("/validate/{message_type}")
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


@router.post("/resolve_input")
async def resolve_input(input_string: str):
    """Resolves an input string into a ReactionInput message."""
    try:
        return send_message(resolvers.resolve_input(input_string))
    except (ValueError, KeyError) as error:
        return Response(str(error), status_code=400)


@router.post("/resolve_compound")
async def resolve_compound(identifier_type: str, compound_string: str):
    """Resolves a compound identifier into a SMILES string."""
    try:
        smiles, resolver = resolvers.name_resolve(identifier_type, compound_string)
        return {"smiles": resolvers.canonicalize_smiles(smiles), "resolver": resolver}
    except ValueError as error:
        return Response(str(error), status_code=400)


@router.post("/canonicalize_smiles")
async def canonicalize_smiles(smiles: str):
    """Canonicalizes a SMILES string."""
    try:
        return resolvers.canonicalize_smiles(smiles)
    except ValueError as error:
        return Response(str(error), status_code=400)


@router.post("/enumerate_dataset")
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
