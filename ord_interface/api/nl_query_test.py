# Copyright 2026 Open Reaction Database Project Authors
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

"""Tests for ord_interface.api.nl_query.

These exercise the serving layer: mapping a translated query onto QueryParams, the
Redis translation cache, and the endpoint. Translation and structure resolution are
tested in ord_schema.agent. The Anthropic client and the name resolver are stubbed, so
no network, API key, or database is required.
"""

import json
from unittest import mock

import pytest
from fastapi import HTTPException
from ord_schema.agent import nl_query as agent_nl_query
from ord_schema.agent.nl_query import NLComponent, NLQuery

from ord_interface.api import nl_query
from ord_interface.api.nl_query import build_query_params
from ord_interface.api.nl_query import (
    nl_query as nl_query_endpoint,
)


@pytest.fixture(autouse=True)
def _no_redis(monkeypatch):
    """Disables Redis so unit tests never touch a real server (always a cache miss)."""

    async def miss(key):
        return None

    async def noop(key, value, ttl_seconds):
        return None

    monkeypatch.setattr(nl_query, "_redis_get", miss)
    monkeypatch.setattr(nl_query, "_redis_set", noop)


@pytest.mark.asyncio
async def test_build_query_params_resolves_name(monkeypatch):
    monkeypatch.setattr(
        agent_nl_query,
        "resolve_name",
        lambda value_type, value: ("CC(=O)O", "PubChem API"),
    )
    query = NLQuery(
        components=[
            NLComponent(identifier="acetic acid", target="INPUT", mode="EXACT")
        ],
        min_yield=70,
    )
    params, resolved = await build_query_params(query)
    assert params.component is not None
    assert json.loads(params.component[0]) == {
        "pattern": "CC(=O)O",
        "target": "INPUT",
        "mode": "EXACT",
    }
    assert params.min_yield == 70
    assert resolved[0].smiles == "CC(=O)O"
    assert resolved[0].resolver == "PubChem API"


@pytest.mark.asyncio
async def test_build_query_params_accepts_verbatim_smiles(monkeypatch):
    # A valid SMILES should never reach the (network) name resolver.
    def fail(*args, **kwargs):
        raise AssertionError("resolver should not be called for a valid SMILES")

    monkeypatch.setattr(agent_nl_query, "resolve_name", fail)
    query = NLQuery(
        components=[
            NLComponent(identifier="c1ccccc1", target="OUTPUT", mode="SUBSTRUCTURE")
        ]
    )
    params, resolved = await build_query_params(query)
    assert params.component is not None
    assert json.loads(params.component[0]) == {
        "pattern": "c1ccccc1",
        "target": "OUTPUT",
        "mode": "SUBSTRUCTURE",
    }
    assert resolved[0].resolver == "SMILES (verbatim)"


@pytest.mark.asyncio
async def test_build_query_params_passes_smarts_through(monkeypatch):
    monkeypatch.setattr(
        agent_nl_query, "canonicalize_smiles", mock.Mock(side_effect=AssertionError)
    )
    query = NLQuery(
        components=[
            NLComponent(identifier="[#6][F,Cl,Br,I]", target="INPUT", mode="SMARTS")
        ]
    )
    params, _ = await build_query_params(query)
    assert params.component is not None
    assert json.loads(params.component[0]) == {
        "pattern": "[#6][F,Cl,Br,I]",
        "target": "INPUT",
        "mode": "SMARTS",
    }


@pytest.mark.asyncio
async def test_build_query_params_invalid_smarts(monkeypatch):
    # A SMARTS the model authored but RDKit cannot parse should be a clean 422, not a
    # 400 surfaced deep in query execution (which a dry run would skip entirely).
    monkeypatch.setattr(
        agent_nl_query, "canonicalize_smiles", mock.Mock(side_effect=AssertionError)
    )
    query = NLQuery(
        components=[NLComponent(identifier="[Br", target="OUTPUT", mode="SMARTS")]
    )
    with pytest.raises(HTTPException) as excinfo:
        await build_query_params(query)
    assert excinfo.value.status_code == 422
    assert "[Br" in excinfo.value.detail


@pytest.mark.asyncio
async def test_build_query_params_invalid_reaction_smarts():
    # An unparseable reaction SMARTS is a clear 422 (in both normal and dry-run mode),
    # not a misleading "no constraints" error from run_query or a silent pass-through.
    query = NLQuery(reaction_smarts="this is not a reaction")
    with pytest.raises(HTTPException) as excinfo:
        await build_query_params(query)
    assert excinfo.value.status_code == 422
    assert "reaction smarts" in excinfo.value.detail.lower()


@pytest.mark.asyncio
async def test_build_query_params_accepts_valid_reaction_smarts():
    query = NLQuery(reaction_smarts="[C:1]=[O:2]>>[C:1][O:2]")
    params, _ = await build_query_params(query)
    assert params.reaction_smarts == "[C:1]=[O:2]>>[C:1][O:2]"


@pytest.mark.asyncio
async def test_build_query_params_unresolvable_name(monkeypatch):
    def raise_value_error(value_type, value):
        raise ValueError(f"Could not resolve {value_type} {value} to SMILES")

    monkeypatch.setattr(agent_nl_query, "resolve_name", raise_value_error)
    query = NLQuery(
        components=[
            NLComponent(identifier="not-a-compound", target="INPUT", mode="EXACT")
        ]
    )
    with pytest.raises(HTTPException) as excinfo:
        await build_query_params(query)
    assert excinfo.value.status_code == 422


@pytest.mark.asyncio
async def test_translation_cache_get_discards_invalid_payload(monkeypatch):
    # A cached entry from an older schema (here, a bad mode) degrades to a miss, not a 500.
    async def stale(key):
        return json.dumps(
            {"components": [{"identifier": "x", "target": "INPUT", "mode": "BOGUS"}]}
        )

    monkeypatch.setattr(nl_query, "_redis_get", stale)
    assert await nl_query._translation_cache_get("key") is None


def _benzene_interpretation() -> NLQuery:
    return NLQuery(
        components=[NLComponent(identifier="benzene", target="INPUT", mode="EXACT")]
    )


@pytest.mark.asyncio
async def test_nl_query_uses_cached_translation_without_model_call(monkeypatch):
    async def fake_cache_get(key):
        return _benzene_interpretation()

    def fail_get_client():
        raise AssertionError("model must not be called on a translation cache hit")

    async def fake_run_query(params, return_ids):
        return []

    monkeypatch.setattr(nl_query, "_translation_cache_get", fake_cache_get)
    monkeypatch.setattr(nl_query, "get_client", fail_get_client)
    monkeypatch.setattr(nl_query, "run_query", fake_run_query)
    monkeypatch.setattr(
        agent_nl_query,
        "resolve_name",
        lambda value_type, value: ("c1ccccc1", "PubChem API"),
    )
    result = await nl_query_endpoint(q="reactions using benzene")
    # The search still runs on a cache hit, so results are fresh.
    assert result.interpretation.components[0].identifier == "benzene"
    assert result.resolved_components[0].smiles == "c1ccccc1"


@pytest.mark.asyncio
async def test_nl_query_translates_and_caches_translation_on_miss(monkeypatch):
    stored = {}

    async def fake_cache_get(key):
        return None

    async def fake_cache_set(key, interpretation):
        stored[key] = interpretation

    async def fake_translate(query, client):
        return _benzene_interpretation()

    async def fake_run_query(params, return_ids):
        return []

    monkeypatch.setattr(nl_query, "_translation_cache_get", fake_cache_get)
    monkeypatch.setattr(nl_query, "_translation_cache_set", fake_cache_set)
    monkeypatch.setattr(nl_query, "get_client", lambda: mock.AsyncMock())
    monkeypatch.setattr(nl_query, "translate", fake_translate)
    monkeypatch.setattr(nl_query, "run_query", fake_run_query)
    monkeypatch.setattr(
        agent_nl_query,
        "resolve_name",
        lambda value_type, value: ("c1ccccc1", "PubChem API"),
    )
    result = await nl_query_endpoint(q="reactions using benzene")
    assert result.resolved_components[0].smiles == "c1ccccc1"
    # Only the translation is cached -- not the search results.
    assert len(stored) == 1
    assert isinstance(next(iter(stored.values())), NLQuery)


@pytest.mark.asyncio
async def test_nl_query_empty_interpretation_returns_422(monkeypatch):
    async def fake_cache_get(key):
        return NLQuery()  # No components, no filters.

    async def fail_run_query(params, return_ids):
        raise ValueError("No query parameters were specified.")

    monkeypatch.setattr(nl_query, "_translation_cache_get", fake_cache_get)
    monkeypatch.setattr(nl_query, "run_query", fail_run_query)
    with pytest.raises(HTTPException) as excinfo:
        await nl_query_endpoint(q="show me everything")
    assert excinfo.value.status_code == 422


@pytest.mark.asyncio
async def test_nl_query_dry_run_skips_search(monkeypatch):
    async def fake_cache_get(key):
        return _benzene_interpretation()

    async def fail_run_query(params, return_ids):
        raise AssertionError("run_query must not be called in dry-run mode")

    monkeypatch.setattr(nl_query, "_translation_cache_get", fake_cache_get)
    monkeypatch.setattr(nl_query, "run_query", fail_run_query)
    monkeypatch.setattr(
        agent_nl_query,
        "resolve_name",
        lambda value_type, value: ("c1ccccc1", "PubChem API"),
    )
    result = await nl_query_endpoint(q="reactions using benzene", dry_run=True)
    assert result.dry_run is True
    assert result.results == []
    # The query that would have run is still surfaced for inspection.
    assert json.loads(result.query_components[0])["pattern"] == "c1ccccc1"


@pytest.mark.parametrize(
    ("error", "expected_status"),
    [
        (nl_query.ModelRateLimitedError("busy"), 429),
        (nl_query.ModelUnavailableError("down"), 503),
        (nl_query.MalformedQueryError("garbled"), 502),
    ],
)
def test_translation_failures_keep_their_status_codes(error, expected_status):
    # ord-schema raises plain exceptions; this mapping is the only thing preserving the
    # API's contract, so it is asserted here rather than inferred from the library.
    http_error = nl_query._as_http_error(error)
    assert http_error.status_code == expected_status
    assert http_error.detail == str(error)


def test_an_unmapped_translation_failure_is_a_500():
    class NewFailure(agent_nl_query.NLQueryError):
        """A subclass added upstream without a status code here."""

    assert nl_query._as_http_error(NewFailure("?")).status_code == 500


@pytest.mark.asyncio
async def test_the_endpoint_surfaces_a_translation_failure(monkeypatch):
    async def fake_cache_get(key):
        return None

    async def failing_translate(query, client):
        raise nl_query.ModelUnavailableError("temporarily unavailable")

    monkeypatch.setattr(nl_query, "_translation_cache_get", fake_cache_get)
    monkeypatch.setattr(nl_query, "get_client", lambda: mock.Mock())
    monkeypatch.setattr(nl_query, "translate", failing_translate)
    with pytest.raises(HTTPException) as excinfo:
        await nl_query_endpoint(q="anything")
    assert excinfo.value.status_code == 503


@pytest.mark.asyncio
async def test_a_missing_api_key_surfaces_as_unavailable(monkeypatch):
    async def fake_cache_get(key):
        return None

    monkeypatch.setattr(nl_query, "_translation_cache_get", fake_cache_get)
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)
    with pytest.raises(HTTPException) as excinfo:
        await nl_query_endpoint(q="anything")
    # get_client() is inside the mapped block, so a missing key is a 503 rather than
    # an unhandled error.
    assert excinfo.value.status_code == 503
