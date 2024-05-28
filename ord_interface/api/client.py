"""Client API."""
import os
from typing import Annotated

from fastapi import APIRouter, Query

from ord_interface.client.query import (
    DatasetIdQuery,
    DoiQuery,
    OrdPostgres,
    ReactionComponentPredicate,
    ReactionComponentQuery,
    ReactionConversionQuery,
    ReactionIdQuery,
    ReactionQuery,
    ReactionSmartsQuery,
    ReactionYieldQuery,
    Result,
)

router = APIRouter(prefix="/api", tags=["client"])

POSTGRES_HOST = os.getenv("POSTGRES_HOST", "localhost")
POSTGRES_PORT = os.getenv("POSTGRES_PORT", "5432")
POSTGRES_USER = os.getenv("POSTGRES_USER", "ord-postgres")
POSTGRES_PASSWORD = os.getenv("POSTGRES_PASSWORD", "ord-postgres")
POSTGRES_DATABASE = os.getenv("POSTGRES_DATABASE", "ord")

BOND_LENGTH = 20
MAX_RESULTS = 1000


def connect():
    return OrdPostgres(
        dbname=POSTGRES_DATABASE,
        user=POSTGRES_USER,
        password=POSTGRES_PASSWORD,
        host=POSTGRES_HOST,
        port=int(POSTGRES_PORT),
    )


def run_query(commands: list[ReactionQuery], limit: int | None) -> list[Result]:
    """Runs a query and returns the matched reactions."""
    if not commands:
        results = []
    elif len(commands) == 1:
        results = connect().run_query(commands[0], limit=limit)
    else:
        # Perform each query without limits and take the intersection of the matched reactions.
        connection = connect()
        reactions = {}
        intersection = None
        for command in commands:
            this_results = {result.reaction_id: result for result in connection.run_query(command)}
            reactions |= this_results
            if intersection is None:
                intersection = set(this_results.keys())
            else:
                intersection &= this_results.keys()
        results = [reactions[reaction_id] for reaction_id in intersection]
        if limit:
            results = results[:limit]
    return results


@router.get("/query")
def query(
    dataset_id: Annotated[list[str] | None, Query()] = None,
    reaction_id: Annotated[list[str] | None, Query()] = None,
    reaction_smarts: str | None = None,
    min_conversion: float | None = None,
    max_conversion: float | None = None,
    min_yield: float | None = None,
    max_yield: float | None = None,
    doi: Annotated[list[str] | None, Query()] = None,
    component: Annotated[list[str] | None, Query()] = None,
    use_stereochemistry: bool | None = None,
    similarity: float | None = None,
    limit: int | None = None,
):
    """Returns a serialized Dataset proto containing matched reactions."""
    queries = []
    if dataset_id:
        queries.append(DatasetIdQuery(dataset_id))
    if reaction_id:
        queries.append(ReactionIdQuery(reaction_id))
    if reaction_smarts:
        queries.append(ReactionSmartsQuery(reaction_smarts))
    if min_conversion is not None or max_conversion is not None:
        queries.append(ReactionConversionQuery(min_conversion, max_conversion))
    if min_yield is not None or max_yield is not None:
        queries.append(ReactionYieldQuery(min_yield, max_yield))
    if doi:
        queries.append(DoiQuery(doi))
    if component:
        predicates = []
        for spec in component:
            pattern, target_name, mode_name = spec.split(";")
            target = ReactionComponentPredicate.Target.from_name(target_name)
            mode = ReactionComponentPredicate.MatchMode.from_name(mode_name)
            predicates.append(ReactionComponentPredicate(pattern, target, mode))
        kwargs = {}
        if use_stereochemistry is not None:
            kwargs["do_chiral_sss"] = use_stereochemistry
        if similarity is not None:
            kwargs["tanimoto_threshold"] = similarity
        queries.append(ReactionComponentQuery(predicates, **kwargs))
    if limit:
        limit = min(limit, MAX_RESULTS)
    else:
        limit = MAX_RESULTS
    if not queries:
        raise
    results = run_query(queries, limit)
    return results
