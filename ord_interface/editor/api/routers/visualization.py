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

"""Visualization API endpoints."""

from fastapi import APIRouter, Request, Response
from ord_schema.message_helpers import mol_from_compound
from ord_schema.proto.reaction_pb2 import Compound, Reaction

from ord_interface.visualization.drawing import mol_to_svg
from ord_interface.visualization.generate_text import generate_html

router = APIRouter(tags=["visualization"])


@router.post("/render_reaction")
async def render_reaction(request: Request):
    """Returns a block of HTML that contains a visual summary of the reaction."""
    reaction = Reaction.FromString(await request.body())
    if not (reaction.inputs or reaction.outcomes):
        return ""
    try:
        return generate_html(reaction)
    except (ValueError, KeyError) as error:
        return Response(str(error), status_code=400)


@router.post("/render_compound")
async def render_compound(request: Request):
    """Returns an HTML-tagged SVG for the given Compound."""
    compound = Compound.FromString(await request.body())
    try:
        return mol_to_svg(mol_from_compound(compound))
    except ValueError as error:
        return Response(str(error), status_code=400)
