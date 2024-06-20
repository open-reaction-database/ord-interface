"""Tests for ord_interface.api.view."""

from ord_schema.proto import reaction_pb2


def test_get_compound_svg(test_client):
    compound = reaction_pb2.Compound()
    compound.identifiers.add(value="c1ccccc1", type="SMILES")
    response = test_client.post("/api/compound_svg", data=compound.SerializeToString())
    response.raise_for_status()


def test_get_reaction_summary(test_client):
    response = test_client.get("/api/reaction_summary", params={"reaction_id": "ord-3f67aa5592fd434d97a577988d3fd241"})
    response.raise_for_status()
