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

from ord_schema import message_helpers, validations
from ord_schema.logging import get_logger
from ord_schema.proto import reaction_pb2
from psycopg import AsyncCursor
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import rdChemReactions

logger = get_logger(__name__)


class ReactionQuery(ABC):
    """Base class for reaction-based queries."""

    @property
    @abstractmethod
    def query_and_parameters(self) -> tuple[str, list]:
        """Returns the query and any query parameters."""


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
            JOIN dataset ON dataset.id = reaction.dataset_id
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
            FROM reaction
            JOIN rdkit.reactions ON rdkit.reactions.id = reaction.rdkit_reaction_id
            WHERE rdkit.reactions.reaction @> reaction_from_smarts(%s)
        """
        return query, [self._reaction_smarts]


class ReactionConversionQuery(ReactionQuery):
    """Looks up reactions by conversion."""

    def __init__(self, min_conversion: float | None, max_conversion: float | None) -> None:
        """Initializes the query.

        Args:
            min_conversion: Minimum conversion, as a percentage.
            max_conversion: Maximum conversion, as a percentage.
        """
        if min_conversion is None and max_conversion is None:
            raise ValueError("At least one of min_conversion or max_conversion must be specified")
        self._min_conversion = min_conversion
        self._max_conversion = max_conversion

    @property
    def query_and_parameters(self) -> tuple[str, list]:
        """Returns the query and any query parameters."""
        query = """
            SELECT DISTINCT reaction.reaction_id
            FROM ord.reaction
            JOIN ord.reaction_outcome on reaction_outcome.reaction_id = reaction.id
            JOIN percentage on percentage.reaction_outcome_id = reaction_outcome.id
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
            JOIN percentage on percentage.product_measurement_id = product_measurement.id
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
            JOIN reaction_provenance ON reaction_provenance.reaction_id = reaction.id
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

    def __init__(  # pylint: disable=too-many-positional-arguments
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
    def query_and_parameters(self) -> tuple[str, list]:
        """Returns the query and any query parameters."""
        if self._target == ReactionComponentQuery.Target.INPUT:
            mols_sql = """
            JOIN reaction_input ON reaction_input.reaction_id = reaction.id
            JOIN compound ON compound.reaction_input_id = reaction_input.id
            JOIN rdkit.mols ON rdkit.mols.id = compound.rdkit_mol_id
            """
        else:
            mols_sql = """
            JOIN reaction_outcome ON reaction_outcome.reaction_id = reaction.id
            JOIN product_compound ON product_compound.reaction_outcome_id = reaction_outcome.id
            JOIN rdkit.mols ON rdkit.mols.id = product_compound.rdkit_mol_id
            """
        if self._match_mode == ReactionComponentQuery.MatchMode.EXACT:
            predicate_sql = "rdkit.mols.smiles = %s"
            params = [Chem.CanonSmiles(self._pattern)]
        elif self._match_mode == ReactionComponentQuery.MatchMode.SIMILAR:
            predicate_sql = "tanimoto_sml(rdkit.mols.morgan_bfp, morganbv_fp(%s)) >= %s"
            params = [self._pattern, self._similarity_threshold]
        elif self._match_mode == ReactionComponentQuery.MatchMode.SUBSTRUCTURE:
            if self._use_chirality:
                predicate_sql = "substruct_chiral(rdkit.mols.mol, %s)"
            else:
                predicate_sql = "substruct(rdkit.mols.mol, %s)"
            params = [self._pattern]
        elif self._match_mode == ReactionComponentQuery.MatchMode.SMARTS:
            if self._use_chirality:
                predicate_sql = "substruct_chiral(rdkit.mols.mol, %s::qmol)"
            else:
                predicate_sql = "substruct(rdkit.mols.mol, %s::qmol)"
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


async def fetch_results(cursor: AsyncCursor) -> list[str]:
    """Fetches query results.

    Args:
        cursor: AsyncCursor instance.

    Returns:
        List of reaction IDs.
    """
    reaction_ids = set()
    async for row in cursor:
        assert row["reaction_id"] not in reaction_ids  # Sanity check for well-written queries.
        reaction_ids.add(row["reaction_id"])
    return list(reaction_ids)


async def run_queries(
    cursor: AsyncCursor,
    reaction_queries: list[ReactionQuery] | ReactionQuery,
    limit: int | None = None,
) -> list[str]:
    """Runs a query against the database.

    Args:
        cursor: AsyncCursor instance.
        reaction_queries: ReactionQuery or list of ReactionQuery.
        limit: Integer maximum number of matches. If None (the default), no limit is set.

    Returns:
        List of reaction IDs.
    """
    if not isinstance(reaction_queries, list):
        reaction_queries = [reaction_queries]
    queries, combined_params = [], []
    for reaction_query in reaction_queries:
        query, params = reaction_query.query_and_parameters
        queries.append(query)
        combined_params.extend(params)
    combined_query = "\nINTERSECT\n".join(queries)
    if limit:
        # TODO(skearnes): Adding LIMIT can significantly slow down queries, especially if they return few results.
        # See https://stackoverflow.com/q/21385555.
        combined_query += "LIMIT %s"
        combined_params.append(limit)
    logger.debug((combined_query, combined_params))
    await cursor.execute(combined_query, combined_params)
    return await fetch_results(cursor)


class QueryResult(BaseModel):
    """Container for a single result from database query."""

    dataset_id: str
    reaction_id: str
    proto: str  # Serialized Reaction protocol buffer (base64).

    @property
    def reaction(self) -> reaction_pb2.Reaction:
        return reaction_pb2.Reaction.FromString(b64decode(self.proto))

    def __eq__(self, other: QueryResult) -> bool:
        return self.dataset_id == other.dataset_id and self.reaction_id == other.reaction_id


async def fetch_reactions(cursor: AsyncCursor, reaction_ids: list[str]) -> list[QueryResult]:
    """Fetches dataset and proto information for a list of reaction IDs."""
    query = """
        SELECT dataset.dataset_id, reaction.reaction_id, reaction.proto
        FROM ord.reaction
        JOIN ord.dataset ON dataset.id = reaction.dataset_id
        WHERE reaction.reaction_id = ANY (%s)
    """
    await cursor.execute(query, (list(set(reaction_ids)),))
    results = []
    async for row in cursor:
        results.append(
            QueryResult(
                dataset_id=row["dataset_id"], reaction_id=row["reaction_id"], proto=b64encode(row["proto"]).decode()
            )
        )
    return results
