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

"""Natural-language query endpoint.

Translation lives in :mod:`ord_schema.agent.nl_query`, next to the schema it describes.
This module is the serving half: it supplies the Redis cache that translation treats as
optional, maps library exceptions onto HTTP status codes, and turns a resolved
interpretation into the :class:`QueryParams` this backend executes.
"""

from __future__ import annotations

import asyncio
import json
from typing import cast

from fastapi import APIRouter, HTTPException, Query
from ord_schema.agent.nl_query import (
    TRANSLATION_CACHE_TTL_SECONDS,
    MalformedQueryError,
    ModelRateLimitedError,
    ModelUnavailableError,
    NLQuery,
    ResolvedComponent,
    get_client,
    resolve_component,
    translate,
    translation_cache_key,
)
from ord_schema.logging import get_logger
from pydantic import BaseModel, ValidationError
from rdkit.Chem import rdChemReactions

from ord_interface.api.queries import QueryResult
from ord_interface.api.search import ComponentSpec, QueryParams, get_redis, run_query

logger = get_logger(__name__)
router = APIRouter(tags=["nl"])

# The cache is an optimization, never a dependency: an unreachable or slow Redis must
# fail fast so the request falls back to a live call instead of stalling on it.
REDIS_OP_TIMEOUT_SECONDS = 1.0

# Model failures are the caller's problem to describe, not the library's. ord-schema
# raises plain exceptions; this is where they become status codes.
_STATUS_CODES = {
    ModelRateLimitedError: 429,
    ModelUnavailableError: 503,
    MalformedQueryError: 502,
}


def _as_http_error(error: Exception) -> HTTPException:
    """Maps a translation failure onto its HTTP status code."""
    for error_type, status_code in _STATUS_CODES.items():
        if isinstance(error, error_type):
            return HTTPException(status_code=status_code, detail=str(error))
    return HTTPException(status_code=500, detail=str(error))


async def _redis_get(key: str) -> str | None:
    """Returns a cached string value, or None on a miss or any Redis failure.

    The cache is best-effort: an unreachable or slow Redis degrades to a miss within
    REDIS_OP_TIMEOUT_SECONDS rather than stalling (or failing) the request.
    """
    try:
        async with asyncio.timeout(REDIS_OP_TIMEOUT_SECONDS):
            async with get_redis() as client:
                value = await client.get(key)
    except Exception as error:
        logger.warning(f"Redis read failed for {key!r}: {error}")
        return None
    if value is None:
        return None
    return value.decode() if isinstance(value, bytes) else value


async def _redis_set(key: str, value: str, ttl_seconds: int) -> None:
    """Stores a string value with a TTL, ignoring an unreachable/slow Redis."""
    try:
        async with asyncio.timeout(REDIS_OP_TIMEOUT_SECONDS):
            async with get_redis() as client:
                await client.set(key, value, ex=ttl_seconds)
    except Exception as error:
        logger.warning(f"Redis write failed for {key!r}: {error}")


class RedisCache:
    """Adapts Redis to ord_schema.agent.nl_query.Cache.

    Both operations are best-effort by contract, which the module-level helpers already
    guarantee: reads degrade to a miss and writes are dropped rather than raising.
    """

    async def get(self, key: str) -> str | None:
        """Returns the cached value for ``key``, or None on a miss or any failure."""
        return await _redis_get(key)

    async def set(self, key: str, value: str, ttl_seconds: int) -> None:
        """Stores ``value`` under ``key``, ignoring failures."""
        await _redis_set(key, value, ttl_seconds)


_CACHE = RedisCache()


async def build_query_params(
    nl_query: NLQuery,
) -> tuple[QueryParams, list[ResolvedComponent]]:
    """Builds structured QueryParams from a translated NLQuery.

    Args:
        nl_query: The model's structured query.

    Returns:
        A tuple of the QueryParams to execute and the resolved components (for display).

    Raises:
        HTTPException: If a named compound cannot be resolved, or a model-authored SMARTS
            or reaction SMARTS is unparseable (all 422).
    """
    if nl_query.reaction_smarts is not None:
        # Validate up front so a bad pattern is a 422 in both normal and dry-run mode
        # (dry run skips run_query, which validates otherwise). ReactionFromSmarts
        # returns None for some malformed input and raises ValueError for the rest.
        try:
            reaction = rdChemReactions.ReactionFromSmarts(nl_query.reaction_smarts)
        except ValueError:
            reaction = None
        if reaction is None:
            raise HTTPException(
                status_code=422,
                detail=f"Invalid reaction SMARTS: {nl_query.reaction_smarts!r}",
            )
    try:
        resolved = await asyncio.gather(
            *(resolve_component(component, _CACHE) for component in nl_query.components)
        )
    except ValueError as error:
        raise HTTPException(status_code=422, detail=str(error)) from error
    components = [
        ComponentSpec(pattern=r.smiles, target=r.target, mode=r.mode).model_dump_json()
        for r in resolved
    ]
    params = QueryParams(
        component=components or None,
        min_yield=nl_query.min_yield,
        max_yield=nl_query.max_yield,
        min_conversion=nl_query.min_conversion,
        max_conversion=nl_query.max_conversion,
        reaction_smarts=nl_query.reaction_smarts,
        similarity=nl_query.similarity_threshold,
        use_stereochemistry=nl_query.use_stereochemistry,
        limit=nl_query.limit,
    )
    return params, resolved


class NLQueryResponse(BaseModel):
    """Response for a natural-language query."""

    query: str
    interpretation: NLQuery
    resolved_components: list[ResolvedComponent]
    # The structured component predicates (JSON ComponentSpec strings) that would be
    # executed -- surfaced so a dry run shows the exact query, not just the results.
    query_components: list[str]
    results: list[QueryResult]
    dry_run: bool = False


async def _translation_cache_get(key: str) -> NLQuery | None:
    """Returns a cached translation, or None on a miss or unparseable payload."""
    raw = await _redis_get(key)
    if raw is None:
        return None
    try:
        return NLQuery.model_validate_json(raw)
    except (ValueError, ValidationError) as error:
        # ValidationError covers a cache entry written by an older NLQuery schema;
        # ValueError covers malformed JSON. Either way, degrade to a miss.
        logger.warning(f"Discarding unparseable cached translation: {error}")
        return None


async def _translation_cache_set(key: str, interpretation: NLQuery) -> None:
    """Stores a translation in the cache (best-effort)."""
    await _redis_set(
        key, interpretation.model_dump_json(), TRANSLATION_CACHE_TTL_SECONDS
    )


@router.get("/nl_query")
async def nl_query(
    q: str = Query(min_length=1, max_length=2000),
    dry_run: bool = False,
) -> NLQueryResponse:
    """Runs a natural-language search.

    The interpreted query and resolved structures are returned alongside the results
    so the user can see -- and trust or correct -- how their question was understood.
    Only the model's translation is cached (best-effort, in Redis): identical questions
    skip the model call, but the database query is always re-run so results stay fresh.
    A Redis outage falls back to a live translation rather than failing the request.

    With ``dry_run=true`` the question is translated and resolved but the database
    search is not executed -- useful for inspecting exactly what query would run.
    """
    key = translation_cache_key(q)
    interpretation = await _translation_cache_get(key)
    if interpretation is not None:
        logger.info(f"NL query translation cache hit for {q!r}")
    else:
        try:
            client = get_client()
            interpretation = await translate(q, client)
        except (ModelUnavailableError, ModelRateLimitedError, MalformedQueryError) as e:
            raise _as_http_error(e) from e
        await _translation_cache_set(key, interpretation)
    # Resolution caches hits and runs the blocking name-resolver lookups in a thread.
    params, resolved = await build_query_params(interpretation)
    logger.info(
        f"NL query {q!r} (dry_run={dry_run}) -> {json.dumps(interpretation.model_dump())}"
    )
    results: list[QueryResult] = []
    if not dry_run:
        try:
            results = cast(list[QueryResult], await run_query(params, return_ids=False))
        except ValueError as error:
            # Raised when the interpretation has no usable constraints (e.g. "show me
            # everything"); surface a 422 rather than an unhandled 500.
            raise HTTPException(
                status_code=422,
                detail="Could not extract any search constraints from the question.",
            ) from error
    return NLQueryResponse(
        query=q,
        interpretation=interpretation,
        resolved_components=resolved,
        query_components=params.component or [],
        results=results,
        dry_run=dry_run,
    )
