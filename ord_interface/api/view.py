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

"""View API."""

from base64 import b64decode

from fastapi import APIRouter
from ord_schema.proto import reaction_pb2
from pydantic import BaseModel

from ord_interface.api.client import ReactionIdList, get_reactions
from ord_interface.visualization import filters, generate_text

router = APIRouter(prefix="/api/view", tags=["view"])


class ViewCompoundInputs(BaseModel):
    """View compound inputs."""

    proto: str  # Serialized Compound protocol buffer (base64).


@router.post("/compound")
def get_compound(inputs: ViewCompoundInputs) -> str:
    """Returns an SVG of a molecule."""
    compound = reaction_pb2.Compound.FromString(b64decode(inputs.proto))
    try:
        return filters._compound_svg(compound)  # pylint: disable=protected-access
    except (ValueError, KeyError):
        return "[Compound cannot be displayed]"


@router.get("/reaction_summary")
def get_reaction_summary(reaction_id: str, compact: bool = True) -> str:
    """Renders a reaction as an HTML table with images and text."""
    results = get_reactions(ReactionIdList(reaction_ids=[reaction_id]))
    if len(results) == 0 or len(results) > 1:
        raise ValueError(reaction_id)
    try:
        return generate_text.generate_html(reaction=results[0].reaction, compact=compact)
    except (ValueError, KeyError):
        return "[Reaction cannot be displayed]"


@router.get("/reaction_page")
def get_reaction_page(reaction_id: str, compact: bool = True):
    """Returns the pbtxt of a single reaction as plain text."""
    results = get_reactions(ReactionIdList(reaction_ids=[reaction_id]))
    if len(results) == 0 or len(results) > 1:
        raise ValueError(reaction_id)
    reaction = results[0].reaction
    reaction_summary = generate_text.generate_html(reaction=reaction, compact=compact)
    return flask.render_template(
        "reaction_view.html",
        reaction=reaction,
        dataset_id=results[0].dataset_id,
        reaction_summary=reaction_summary,
        bond_length=20,
    )
