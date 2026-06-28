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

"""Natural-language query interface.

Translates a free-text question (e.g. "find reactions for synthesizing ibuprofen
with yield greater than 70%") into the structured ``QueryParams`` understood by
the existing search backend, then dispatches it through the same index-accelerated
path as the structured API.

The language model only ever emits a structured query (a forced tool call); it
never writes SQL and never invents SMILES. Compound names are resolved to SMILES
deterministically via ``ord_schema.resolvers``, so the model's chemistry is
grounded in PubChem/CIR/OPSIN rather than its own recall.
"""

from __future__ import annotations

import asyncio
import hashlib
import json
import os
from importlib import resources
from typing import Literal, cast

import anthropic
from fastapi import APIRouter, HTTPException, Query
from ord_schema.logging import get_logger
from ord_schema.resolvers import canonicalize_smiles, resolve_name
from pydantic import BaseModel, Field, ValidationError
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from ord_interface.api.queries import QueryResult
from ord_interface.api.search import ComponentSpec, QueryParams, get_redis, run_query

logger = get_logger(__name__)
router = APIRouter(tags=["nl"])

# Haiku is fast and cheap and the translation task is tightly constrained; override
# with ORD_NL_QUERY_MODEL (e.g. a Sonnet snapshot) if disambiguation needs more reasoning.
DEFAULT_MODEL = "claude-haiku-4-5"
MAX_TOKENS = 1024

# Only the model's translation is cached -- not the search results -- so repeated
# identical questions don't each pay for a model call, while the database query is always
# re-run and stays fresh. Bump the version when the prompt or NLQuery schema changes so
# stale interpretations are not served.
TRANSLATION_CACHE_VERSION = "v3"
TRANSLATION_CACHE_TTL_SECONDS = 60 * 60

# Name -> SMILES resolutions are cached separately and for much longer: they are stable
# (a name maps to the same structure) and shared across different questions that mention
# the same compound, which spares the name resolvers repeated lookups.
RESOLVE_CACHE_VERSION = "v1"
RESOLVE_CACHE_TTL_SECONDS = 60 * 60 * 24 * 30

# The cache is an optimization, never a dependency: an unreachable or slow Redis must
# fail fast so the request falls back to a live call instead of stalling on it.
REDIS_OP_TIMEOUT_SECONDS = 1.0

Target = Literal["INPUT", "OUTPUT"]
MatchMode = Literal["EXACT", "SIMILAR", "SUBSTRUCTURE", "SMARTS"]


class NLComponent(BaseModel):
    """A single component constraint extracted from the question."""

    identifier: str = Field(
        description=(
            "The compound as the user named it: a common/trade/IUPAC name (e.g. "
            "'ibuprofen', 'benzene'), a SMILES, or -- only when mode is SMARTS -- a "
            "SMARTS pattern. Prefer the plain name; do not translate names to SMILES."
        )
    )
    target: Target = Field(
        description=(
            "INPUT if the compound is a reactant/reagent/solvent consumed by the "
            "reaction ('using X', 'from X', 'with X'); OUTPUT if it is produced by "
            "the reaction ('synthesizing X', 'to make X', 'yields X')."
        )
    )
    mode: MatchMode = Field(
        description=(
            "EXACT for a specific named molecule (the default). SUBSTRUCTURE when the "
            "user wants molecules that merely contain a group/scaffold ('containing a "
            "pyridine ring'). SIMILAR for 'similar to'/'like X'. SMARTS only when the "
            "user supplies or describes an explicit query pattern."
        )
    )


class NLQuery(BaseModel):
    """Structured query produced by the language model from natural language."""

    components: list[NLComponent] = Field(
        default_factory=list,
        description="Per-compound constraints; AND-combined with the other fields.",
    )
    min_yield: float | None = Field(
        default=None, description="Minimum percent yield (0-100), if requested."
    )
    max_yield: float | None = Field(
        default=None, description="Maximum percent yield (0-100), if requested."
    )
    min_conversion: float | None = Field(
        default=None, description="Minimum percent conversion (0-100), if requested."
    )
    max_conversion: float | None = Field(
        default=None, description="Maximum percent conversion (0-100), if requested."
    )
    reaction_smarts: str | None = Field(
        default=None,
        description=(
            "A reaction SMARTS (reactants>>products) when the user describes a "
            "transformation rather than individual components. Usually omitted."
        ),
    )
    similarity_threshold: float | None = Field(
        default=None,
        description="Tanimoto threshold (0-1) for SIMILAR components; default 0.5.",
    )
    use_stereochemistry: bool | None = Field(
        default=None,
        description="True only if the user asks to respect stereochemistry/chirality.",
    )
    limit: int | None = Field(
        default=None, description="Maximum number of reactions to return, if stated."
    )


# The system prompt lives in nl_query_prompt.md (alongside this module) so it reads as
# plain markdown and can be edited without touching Python string escaping.
SYSTEM_PROMPT = (
    (resources.files("ord_interface.api") / "nl_query_prompt.md")
    .read_text(encoding="utf-8")
    .strip()
)

_TOOL = {
    "name": "build_query",
    "description": "Build a structured ORD search query from the user's question.",
    "input_schema": NLQuery.model_json_schema(),
}


def _get_client() -> anthropic.AsyncAnthropic:
    """Returns an Anthropic async client.

    Raises:
        HTTPException: If ANTHROPIC_API_KEY is not configured.
    """
    if not os.getenv("ANTHROPIC_API_KEY"):
        raise HTTPException(
            status_code=503,
            detail="Natural-language search is unavailable: ANTHROPIC_API_KEY is not set.",
        )
    return anthropic.AsyncAnthropic()


async def translate(query: str, client: anthropic.AsyncAnthropic) -> NLQuery:
    """Translates a natural-language question into a structured query.

    Args:
        query: The user's free-text question.
        client: Anthropic async client.

    Returns:
        The structured NLQuery produced by the model.

    Raises:
        HTTPException: If the model is rate limited (429), otherwise unreachable or
            erroring (503), or does not return a usable structured query (502).
    """
    try:
        response = await client.messages.create(
            model=os.getenv("ORD_NL_QUERY_MODEL", DEFAULT_MODEL),
            max_tokens=MAX_TOKENS,
            system=SYSTEM_PROMPT,
            tools=[cast(anthropic.types.ToolParam, _TOOL)],
            tool_choice={"type": "tool", "name": "build_query"},
            messages=[{"role": "user", "content": query}],
        )
    except anthropic.RateLimitError as error:
        raise HTTPException(
            status_code=429,
            detail="Natural-language search is busy right now; please retry shortly.",
        ) from error
    except anthropic.APIError as error:
        # Connection failures, server errors, overloads, and auth problems all degrade
        # to a graceful "temporarily unavailable" rather than a 500.
        logger.warning(f"Anthropic API error during NL translation: {error}")
        raise HTTPException(
            status_code=503,
            detail="Natural-language search is temporarily unavailable.",
        ) from error
    for block in response.content:
        if (
            isinstance(block, anthropic.types.ToolUseBlock)
            and block.name == "build_query"
        ):
            try:
                return NLQuery.model_validate(block.input)
            except ValidationError as error:
                # The forced tool schema makes this unlikely, but a payload that slips
                # through becomes a 502 rather than an unhandled 500.
                logger.warning(f"Model tool call failed schema validation: {error}")
                raise HTTPException(
                    status_code=502,
                    detail="Language model returned a malformed structured query.",
                ) from error
    raise HTTPException(
        status_code=502, detail="Language model did not return a structured query."
    )


class ResolvedComponent(BaseModel):
    """A component after name resolution, surfaced for transparency."""

    identifier: str
    smiles: str
    resolver: str
    target: Target
    mode: MatchMode


def _resolve_name_key(name: str) -> str:
    """Returns the Redis cache key for a name -> SMILES resolution."""
    digest = hashlib.sha256(name.strip().lower().encode()).hexdigest()
    return f"nl_resolve:{RESOLVE_CACHE_VERSION}:{digest}"


async def _resolve_name_cached(name: str) -> tuple[str, str]:
    """Resolves a compound name to (SMILES, resolver), caching successful lookups.

    The blocking PubChem/CIR/OPSIN call runs in a worker thread. Failures are not cached so
    a transient PubChem outage does not poison the cache with a permanent miss.

    Args:
        name: The compound name to resolve.

    Returns:
        A tuple of canonical SMILES and the resolver that produced it (e.g. "PubChem
        API"); the resolver is suffixed with " (cached)" on a cache hit.

    Raises:
        ValueError: If the name cannot be resolved to a structure.
    """
    key = _resolve_name_key(name)
    raw = await _redis_get(key)
    if raw is not None:
        try:
            smiles, resolver = json.loads(raw)
            return smiles, f"{resolver} (cached)"
        except (ValueError, TypeError) as error:
            logger.warning(f"Discarding bad cached resolution for {name!r}: {error}")
    smiles, resolver = await asyncio.to_thread(resolve_name, "name", name)
    await _redis_set(key, json.dumps([smiles, resolver]), RESOLVE_CACHE_TTL_SECONDS)
    return smiles, resolver


async def _resolve_component(component: NLComponent) -> ResolvedComponent:
    """Resolves a component's identifier to canonical SMILES.

    SMARTS patterns pass through untouched once validated. Otherwise the identifier is
    treated as a SMILES if RDKit can parse it, and falls back to (cached) name resolution
    via PubChem/CIR/OPSIN.

    Args:
        component: The component to resolve.

    Returns:
        The resolved component.

    Raises:
        ValueError: If a SMARTS pattern is unparseable, or a non-SMARTS identifier can be
            resolved to neither SMILES nor a name.
    """
    if component.mode == "SMARTS":
        # The model authors SMARTS directly; validate up front so a bad pattern is a
        # 422 here rather than a 400 deep in query execution (skipped on dry runs).
        if Chem.MolFromSmarts(component.identifier) is None:
            raise ValueError(f"Invalid SMARTS pattern: {component.identifier!r}")
        return ResolvedComponent(
            identifier=component.identifier,
            smiles=component.identifier,
            resolver="SMARTS (verbatim)",
            target=component.target,
            mode=component.mode,
        )
    try:
        smiles = canonicalize_smiles(component.identifier)
        resolver = "SMILES (verbatim)"
    except ValueError:
        smiles, resolver = await _resolve_name_cached(component.identifier)
    return ResolvedComponent(
        identifier=component.identifier,
        smiles=smiles,
        resolver=resolver,
        target=component.target,
        mode=component.mode,
    )


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
            *(_resolve_component(component) for component in nl_query.components)
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
    """Stores a string value with a TTL, ignoring an unreachable/slow Redis (best-effort)."""
    try:
        async with asyncio.timeout(REDIS_OP_TIMEOUT_SECONDS):
            async with get_redis() as client:
                await client.set(key, value, ex=ttl_seconds)
    except Exception as error:
        logger.warning(f"Redis write failed for {key!r}: {error}")


def _translation_cache_key(query: str) -> str:
    """Returns the Redis cache key for a question under the current model and version."""
    model = os.getenv("ORD_NL_QUERY_MODEL", DEFAULT_MODEL)
    digest = hashlib.sha256(f"{model}\n{query.strip()}".encode()).hexdigest()
    return f"nl_query:{TRANSLATION_CACHE_VERSION}:{digest}"


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
    key = _translation_cache_key(q)
    interpretation = await _translation_cache_get(key)
    if interpretation is not None:
        logger.info(f"NL query translation cache hit for {q!r}")
    else:
        client = _get_client()
        interpretation = await translate(q, client)
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
