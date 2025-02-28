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

"""Tests for ord_interface.api.view."""

from ord_schema.proto import reaction_pb2


def test_get_compound_svg(test_client):
    compound = reaction_pb2.Compound()
    compound.identifiers.add(value="c1ccccc1", type="SMILES")
    response = test_client.post("/api/compound_svg", content=compound.SerializeToString())
    response.raise_for_status()


def test_get_reaction_summary(test_client):
    response = test_client.get("/api/reaction_summary", params={"reaction_id": "ord-3f67aa5592fd434d97a577988d3fd241"})
    response.raise_for_status()
