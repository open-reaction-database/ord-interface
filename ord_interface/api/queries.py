# Copyright 2020 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Library for executing PostgreSQL queries on the ORD.

A reaction query consists of _one_ of the following:

    * A reaction ID
    * A reaction SMARTS
    * A set of reaction component predicates.

Each reaction component predicate has the following structure:

    * Input/output selector
    * One of the following:
        * Exact structure match
        * Substructure match (including SMARTS)
        * Structural similarity

    Note that similarity searches use a query-level similarity threshold; it is
    not possible to set predicate-level thresholds (unless the predicates are
    run as separate queries or some sort of post-hoc filtering is used).

For example, a reaction query might have the following predicates:

    * Input is c1ccccc1 (exact match)
    * Input contains C(=O)O (substructure search)
    * Output has at least 0.6 Tanimoto similarity to O=C(C)Oc1ccccc1C(=O)O

Note that a predicate is matched if it applies to _any_ input/output.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from base64 import b64decode, b64encode
from enum import Enum, auto
from typing import Any, LiteralString

from ord_schema import message_helpers, validations
from ord_schema.logging import get_logger
from ord_schema.proto import reaction_pb2
from psycopg import AsyncCursor, sql
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import rdChemReactions

DictCursor = AsyncCursor[dict[str, Any]]

logger = get_logger(__name__)


class ReactionQuery(ABC):
    """Base class for reaction-based queries."""

    @property
    @abstractmethod
    def query_and_parameters(self) -> tuple[str, list]:
        """Returns the query and any query parameters."""

    @property
    def session_config(self) -> dict[str, str]:
        """RDKit GUC settings required for index-assisted execution.

        Some RDKit operators read tuning parameters from session settings rather
        than from the query text (e.g. the similarity threshold for ``%`` and the
        chirality flag for ``@>``). ``run_queries`` applies these via
        ``set_config`` before executing.

        Returns:
            Mapping of GUC name to value, as strings for ``set_config()``.
        """
        return {}


class DatasetIdQuery(ReactionQuery):
    """Looks up reactions by dataset ID."""

    def __init__(self, dataset_ids: list[str]) -> None:
        """Initializes the query.

        Args:
            dataset_ids: List of dataset IDs.
        """
        for dataset_id in dataset_ids:
            if not validations.is_valid_dataset_id(dataset_id):
                raise ValueError(f"Invalid dataset ID: {dataset_id}")
        self._dataset_ids = dataset_ids

    @property
    def query_and_parameters(self) -> tuple[str, list]:
        """Returns the query and any query parameters."""
        query = """
            SELECT DISTINCT reaction.reaction_id
            FROM ord.reaction
            JOIN ord.dataset ON dataset.id = reaction.dataset_id
            WHERE dataset.dataset_id = ANY (%s)
        """
        return query, [self._dataset_ids]


class ReactionIdQuery(ReactionQuery):
    """Looks up reactions by ID."""

    def __init__(self, reaction_ids: list[str]) -> None:
        """Initializes the query.

        Args:
            reaction_ids: List of reaction IDs.
        """
        for reaction_id in reaction_ids:
            if not validations.is_valid_reaction_id(reaction_id):
                raise ValueError(f"Invalid reaction ID: {reaction_id}")
        self._reaction_ids = reaction_ids

    @property
    def query_and_parameters(self) -> tuple[str, list]:
        """Returns the query and any query parameters."""
        query = """
            SELECT DISTINCT reaction.reaction_id
            FROM ord.reaction
            WHERE reaction.reaction_id = ANY (%s)
        """
        return query, [self._reaction_ids]


class ReactionSmartsQuery(ReactionQuery):
    """Matches reactions by reaction SMARTS."""

    def __init__(self, reaction_smarts: str) -> None:
        """Initializes the query.

        Args:
            reaction_smarts: Reaction SMARTS.
        """
        reaction = rdChemReactions.ReactionFromSmarts(reaction_smarts)
        if not reaction:
            raise ValueError(f"Cannot parse reaction SMARTS: {reaction_smarts}")
        self._reaction_smarts = reaction_smarts

    @property
    def query_and_parameters(self) -> tuple[str, list]:
        """Returns the query and any query parameters."""
        query = """
            SELECT DISTINCT reaction.reaction_id
            FROM ord.reaction
            JOIN rdkit.reactions ON rdkit.reactions.id = reaction.rdkit_reaction_id
            WHERE rdkit.reactions.reaction @> reaction_from_smarts(%s)
        """
        return query, [self._reaction_smarts]


class ReactionConversionQuery(ReactionQuery):
    """Looks up reactions by conversion."""

    def __init__(
        self, min_conversion: float | None, max_conversion: float | None
    ) -> None:
        """Initializes the query.

        Args:
            min_conversion: Minimum conversion, as a percentage.
            max_conversion: Maximum conversion, as a percentage.
        """
        if min_conversion is None and max_conversion is None:
            raise ValueError(
                "At least one of min_conversion or max_conversion must be specified"
            )
        self._min_conversion = min_conversion
        self._max_conversion = max_conversion

    @property
    def query_and_parameters(self) -> tuple[str, list]:
        """Returns the query and any query parameters."""
        query = """
            SELECT DISTINCT reaction.reaction_id
            FROM ord.reaction
            JOIN ord.reaction_outcome on reaction_outcome.reaction_id = reaction.id
            JOIN ord.percentage on percentage.reaction_outcome_id = reaction_outcome.id
        """
        if self._min_conversion is not None and self._max_conversion is not None:
            query += "WHERE percentage.value >= %s AND percentage.value <= %s\n"
            params = [self._min_conversion, self._max_conversion]
        elif self._min_conversion is not None:
            query += "WHERE percentage.value >= %s\n"
            params = [self._min_conversion]
        else:
            query += "WHERE percentage.value <= %s\n"
            params = [self._max_conversion]
        return query, params


class ReactionYieldQuery(ReactionQuery):
    """Looks up reactions by yield."""

    def __init__(self, min_yield: float | None, max_yield: float | None) -> None:
        """Initializes the query.

        Args:
            min_yield: Minimum yield, as a percentage.
            max_yield: Maximum yield, as a percentage.
        """
        if min_yield is None and max_yield is None:
            raise ValueError("At least one of min_yield or max_yield must be specified")
        self._min_yield = min_yield
        self._max_yield = max_yield

    @property
    def query_and_parameters(self) -> tuple[str, list]:
        """Returns the query and any query parameters."""
        query = """
            SELECT DISTINCT reaction.reaction_id
            FROM ord.reaction
            JOIN ord.reaction_outcome on reaction_outcome.reaction_id = reaction.id
            JOIN ord.product_compound on product_compound.reaction_outcome_id = reaction_outcome.id
            JOIN ord.product_measurement on product_measurement.product_compound_id = product_compound.id
            JOIN ord.percentage on percentage.product_measurement_id = product_measurement.id
            WHERE product_measurement.type = 'YIELD'
        """
        params = []
        if self._min_yield is not None:
            query += "AND percentage.value >= %s\n"
            params.append(self._min_yield)
        if self._max_yield is not None:
            query += "AND percentage.value <= %s\n"
            params.append(self._max_yield)
        return query, params


class DoiQuery(ReactionQuery):
    """Looks up reactions by DOI."""

    def __init__(self, dois: list[str]) -> None:
        """Initializes the query.

        Args:
            dois: List of DOIs.
        """
        parsed_dois = []
        for doi in dois:
            try:
                parsed = message_helpers.parse_doi(doi)
            except ValueError as error:
                raise ValueError(f"Invalid DOI: {doi}") from error
            if doi != parsed:
                # Trim the DOI as needed to match the database contents.
                logger.debug(f"Updating DOI: {doi} -> {parsed}")
            parsed_dois.append(parsed)
        self._dois = parsed_dois

    @property
    def query_and_parameters(self) -> tuple[str, list]:
        """Returns the query and any query parameters."""
        query = """
            SELECT DISTINCT reaction.reaction_id
            FROM ord.reaction
            JOIN ord.reaction_provenance ON reaction_provenance.reaction_id = reaction.id
            WHERE reaction_provenance.doi = ANY (%s)
        """
        return query, [self._dois]


class ReactionComponentQuery(ReactionQuery):
    """Matches reactions according to a component-level search."""

    class Target(Enum):
        """Search targets."""

        INPUT = auto()
        OUTPUT = auto()

    class MatchMode(Enum):
        """Interpretations for SMILES and SMARTS strings."""

        EXACT = auto()
        SIMILAR = auto()
        SUBSTRUCTURE = auto()
        SMARTS = auto()

    def __init__(
        self,
        pattern: str,
        target: Target,
        match_mode: MatchMode,
        similarity_threshold: float = 0.5,
        use_chirality: bool = False,
    ) -> None:
        """Initializes the predicate.

        Args:
            pattern: SMILES or SMARTS pattern.
            target: ReactionComponentQuery.Target.
            match_mode: ReactionComponentQuery.MatchMode.
            similarity_threshold: Similarity threshold for SIMILAR mode.
            use_chirality: Whether to use chirality in SUBSTRUCTURE/SMARTS modes.
        """
        if match_mode == ReactionComponentQuery.MatchMode.SMARTS:
            mol = Chem.MolFromSmarts(pattern)
        else:
            mol = Chem.MolFromSmiles(pattern)
        if not mol:
            raise ValueError(f"Cannot parse pattern: {pattern}")
        self._pattern = pattern
        self._target = target
        self._match_mode = match_mode
        self._similarity_threshold = similarity_threshold
        self._use_chirality = use_chirality

    @property
    def _mols_join(self) -> LiteralString:
        """Returns the JOINs linking ord.reaction to rdkit.mols for the target components."""
        if self._target == ReactionComponentQuery.Target.INPUT:
            return """
            JOIN ord.reaction_input ON reaction_input.reaction_id = reaction.id
            JOIN ord.compound ON compound.reaction_input_id = reaction_input.id
            JOIN rdkit.mols ON rdkit.mols.id = compound.rdkit_mol_id
            """
        return """
            JOIN ord.reaction_outcome ON reaction_outcome.reaction_id = reaction.id
            JOIN ord.product_compound ON product_compound.reaction_outcome_id = reaction_outcome.id
            JOIN rdkit.mols ON rdkit.mols.id = product_compound.rdkit_mol_id
            """

    @property
    def is_similarity(self) -> bool:
        """Whether this is a similarity query, whose matches can be ranked by score."""
        return self._match_mode == ReactionComponentQuery.MatchMode.SIMILAR

    @property
    def query_and_parameters(self) -> tuple[str, list]:
        """Returns the query and any query parameters."""
        mols_sql = self._mols_join
        if self._match_mode == ReactionComponentQuery.MatchMode.EXACT:
            predicate_sql = "rdkit.mols.smiles = %s"
            params = [Chem.CanonSmiles(self._pattern)]
        elif self._match_mode == ReactionComponentQuery.MatchMode.SIMILAR:
            # The %% (Tanimoto) operator uses the GiST index on morgan_bfp; the
            # threshold is supplied via rdkit.tanimoto_threshold (session_config).
            predicate_sql = "rdkit.mols.morgan_bfp %% morganbv_fp(%s)"
            params = [self._pattern]
        elif self._match_mode == ReactionComponentQuery.MatchMode.SUBSTRUCTURE:
            # The @> (substructure) operator uses the GiST index on mol; chirality
            # is controlled via rdkit.do_chiral_sss (session_config).
            predicate_sql = "rdkit.mols.mol @> %s::mol"
            params = [self._pattern]
        elif self._match_mode == ReactionComponentQuery.MatchMode.SMARTS:
            predicate_sql = "rdkit.mols.mol @> %s::qmol"
            params = [self._pattern]
        else:
            raise NotImplementedError(f"Unsupported match_mode: {self._match_mode}")
        query = f"""
            SELECT DISTINCT reaction.reaction_id
            FROM ord.reaction
            {mols_sql}
            WHERE {predicate_sql}
        """
        return query, params

    def similarity_score_query(self) -> tuple[LiteralString, list]:
        """Returns a query ranking reactions by similarity, plus its leading parameter.

        Each reaction is scored by the greatest Tanimoto similarity of its target
        components to the query pattern, and the results are ordered by that score,
        descending. The caller supplies the candidate reaction IDs (and any LIMIT) as
        the remaining parameters; scoring runs over that already-filtered set, so the
        per-row ``tanimoto_sml()`` calls are cheap.

        Raises:
            ValueError: If this is not a similarity query.
        """
        if not self.is_similarity:
            raise ValueError("similarity_score_query is only valid for SIMILAR queries")
        query = f"""
            SELECT reaction.reaction_id,
                   MAX(tanimoto_sml(rdkit.mols.morgan_bfp, morganbv_fp(%s))) AS similarity
            FROM ord.reaction
            {self._mols_join}
            WHERE reaction.reaction_id = ANY (%s)
            GROUP BY reaction.reaction_id
            ORDER BY similarity DESC, reaction.reaction_id
        """
        return query, [self._pattern]

    @property
    def session_config(self) -> dict[str, str]:
        """RDKit GUC settings required for index-assisted execution."""
        if self._match_mode == ReactionComponentQuery.MatchMode.SIMILAR:
            return {"rdkit.tanimoto_threshold": str(self._similarity_threshold)}
        if (
            self._match_mode
            in (
                ReactionComponentQuery.MatchMode.SUBSTRUCTURE,
                ReactionComponentQuery.MatchMode.SMARTS,
            )
            and self._use_chirality
        ):
            return {"rdkit.do_chiral_sss": "true"}
        return {}


async def fetch_results(cursor: DictCursor) -> list[str]:
    """Fetches query results.

    Args:
        cursor: AsyncCursor instance.

    Returns:
        List of reaction IDs.
    """
    reaction_ids = set()
    async for row in cursor:
        assert (
            row["reaction_id"] not in reaction_ids
        )  # Sanity check for well-written queries.
        reaction_ids.add(row["reaction_id"])
    return list(reaction_ids)


async def run_queries(
    cursor: DictCursor,
    reaction_queries: list[ReactionQuery] | ReactionQuery,
    limit: int | None = None,
) -> list[str]:
    """Runs a query against the database.

    Args:
        cursor: AsyncCursor instance.
        reaction_queries: ReactionQuery or list of ReactionQuery.
        limit: Integer maximum number of matches. If None (the default), no limit is set.

    Returns:
        List of reaction IDs. When the query contains exactly one similarity
        predicate, the IDs are ordered by descending similarity; otherwise the
        order is unspecified.
    """
    queries_list: list[ReactionQuery] = (
        [reaction_queries]
        if isinstance(reaction_queries, ReactionQuery)
        else reaction_queries
    )
    similarity_queries = [
        query
        for query in queries_list
        if isinstance(query, ReactionComponentQuery) and query.is_similarity
    ]
    # Ranking needs a single similarity predicate to define the order; with zero or
    # several, results are returned unordered.
    ranking_query = similarity_queries[0] if len(similarity_queries) == 1 else None

    queries, combined_params = [], []
    config: dict[str, str] = {}
    for reaction_query in queries_list:
        query, params = reaction_query.query_and_parameters
        queries.append(query)
        combined_params.extend(params)
        for name, value in reaction_query.session_config.items():
            existing = config.setdefault(name, value)
            if existing != value:
                raise ValueError(
                    f"Conflicting values for {name}: {existing} != {value}"
                )
    combined_query = "\nINTERSECT\n".join(queries)
    if limit and ranking_query is None:
        # LIMIT sits on top of the whole INTERSECT. Each branch is materialized
        # (sort + unique for DISTINCT) before the set operation, so this truncates
        # the final result rather than bounding the per-branch index scans. When
        # ranking, the LIMIT is deferred to the scoring step so it selects the most
        # similar matches rather than an arbitrary subset.
        combined_query += "LIMIT %s"
        combined_params.append(limit)
    # Apply RDKit GUCs so the substructure/similarity operators can use their GiST
    # indexes. Connections are created per request (see get_cursor), so these
    # session-level settings do not leak between queries.
    for name, value in config.items():
        await cursor.execute("SELECT set_config(%s, %s, false)", (name, value))
    logger.debug((combined_query, combined_params))
    await cursor.execute(combined_query, combined_params)
    reaction_ids = await fetch_results(cursor)
    if ranking_query is not None:
        reaction_ids = await _rank_by_similarity(
            cursor, ranking_query, reaction_ids, limit
        )
    return reaction_ids


async def _rank_by_similarity(
    cursor: DictCursor,
    ranking_query: ReactionComponentQuery,
    reaction_ids: list[str],
    limit: int | None,
) -> list[str]:
    """Returns the matched reaction IDs ordered by descending similarity, limited to ``limit``."""
    if not reaction_ids:
        return []
    query, params = ranking_query.similarity_score_query()
    params.append(reaction_ids)
    if limit:
        query += "\nLIMIT %s"
        params.append(limit)
    await cursor.execute(query, params)
    return [row["reaction_id"] async for row in cursor]


class QueryResult(BaseModel):
    """Container for a single result from database query."""

    dataset_id: str
    reaction_id: str
    proto: str  # Serialized Reaction protocol buffer (base64).

    @property
    def reaction(self) -> reaction_pb2.Reaction:
        return reaction_pb2.Reaction.FromString(b64decode(self.proto))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, QueryResult):
            return NotImplemented
        return (
            self.dataset_id == other.dataset_id
            and self.reaction_id == other.reaction_id
        )


async def fetch_reactions(
    cursor: DictCursor, reaction_ids: list[str]
) -> list[QueryResult]:
    """Fetches dataset and proto information for a list of reaction IDs.

    Results preserve the order of ``reaction_ids`` (e.g. similarity ranking from
    ``run_queries``); duplicate or unknown IDs are dropped.
    """
    unique_ids = list(dict.fromkeys(reaction_ids))  # De-dup, preserving order.
    query = """
        SELECT dataset.dataset_id, reaction.reaction_id, reaction.proto
        FROM ord.reaction
        JOIN ord.dataset ON dataset.id = reaction.dataset_id
        WHERE reaction.reaction_id = ANY (%s)
    """
    await cursor.execute(query, (unique_ids,))
    results_by_id = {}
    async for row in cursor:
        results_by_id[row["reaction_id"]] = QueryResult(
            dataset_id=row["dataset_id"],
            reaction_id=row["reaction_id"],
            proto=b64encode(row["proto"]).decode(),
        )
    return [results_by_id[rid] for rid in unique_ids if rid in results_by_id]


class StatsResult(BaseModel):
    """Stats result."""

    smiles: str
    times_appearing: int


async def _fetch_dataset_most_used_smiles(
    cursor: DictCursor,
    dataset_id: str,
    *,
    compound_table: str,
    join_table: str,
    foreign_key: str,
    limit: int,
) -> list[StatsResult]:
    """Fetches the top K most used SMILES for a dataset, joined through the given tables.

    Args:
        cursor: Database cursor.
        dataset_id: Dataset to aggregate over.
        compound_table: Table holding SMILES, e.g. "compound" or "product_compound".
        join_table: Table linking compounds to reactions, e.g. "reaction_input".
        foreign_key: Column on ``compound_table`` referencing ``join_table``.
        limit: Maximum number of rows to return.

    Returns:
        Most frequently appearing SMILES, in descending order of frequency.
    """
    # Compose schema-qualified identifiers safely rather than interpolating raw
    # strings into the SQL text.
    compound = sql.Identifier("ord", compound_table)
    join = sql.Identifier("ord", join_table)
    query = sql.SQL(
        """
            SELECT smiles, COUNT(*) AS times_appearing
            FROM {compound}
            JOIN {join} ON {compound}.{foreign_key} = {join}.id
            JOIN ord.reaction ON {join}.reaction_id = ord.reaction.id
            JOIN ord.dataset ON ord.reaction.dataset_id = ord.dataset.id
            WHERE ord.dataset.dataset_id = %s
            AND smiles IS NOT NULL
            GROUP BY smiles
            ORDER BY times_appearing DESC
            LIMIT %s
        """
    ).format(compound=compound, join=join, foreign_key=sql.Identifier(foreign_key))
    await cursor.execute(query, (dataset_id, limit))
    results = []
    async for row in cursor:
        results.append(StatsResult(**row))
    return results


async def fetch_dataset_most_used_smiles_for_inputs(
    cursor: DictCursor, dataset_id: str, limit: int = 30
) -> list[StatsResult]:
    """Fetches the top K most used SMILES molecules in terms of reaction inputs for a given dataset."""
    return await _fetch_dataset_most_used_smiles(
        cursor,
        dataset_id,
        compound_table="compound",
        join_table="reaction_input",
        foreign_key="reaction_input_id",
        limit=limit,
    )


async def fetch_dataset_most_used_smiles_for_products(
    cursor: DictCursor, dataset_id: str, limit: int = 30
) -> list[StatsResult]:
    """Fetches the top K most used SMILES molecules in terms of reaction products for a given dataset."""
    return await _fetch_dataset_most_used_smiles(
        cursor,
        dataset_id,
        compound_table="product_compound",
        join_table="reaction_outcome",
        foreign_key="reaction_outcome_id",
        limit=limit,
    )
