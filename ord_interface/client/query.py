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

import abc
import dataclasses
import enum
import json
import logging
from typing import Dict, List, Optional, Tuple

import psycopg2
import psycopg2.extensions
from psycopg2 import sql
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import reaction_pb2

logger = logging.getLogger()


@dataclasses.dataclass(frozen=True)
class Result:
    """Container for a single result from database query."""

    dataset_id: str
    reaction_id: str
    proto: Optional[bytes] = None

    @property
    def reaction(self) -> reaction_pb2.Reaction:
        return reaction_pb2.Reaction.FromString(self.proto)

    def __eq__(self, other: Result) -> bool:
        return self.dataset_id == other.dataset_id and self.reaction_id == other.reaction_id


def fetch_results(cursor: psycopg2.extensions.cursor) -> List[Result]:
    """Fetches query results.

    Args:
        cursor: psycopg.cursor instance.

    Returns:
        List of Result instances.
    """
    results = []
    for dataset_id, reaction_id, proto in cursor:
        results.append(Result(reaction_id=reaction_id, dataset_id=dataset_id, proto=proto.tobytes()))
    return results


class ReactionQueryBase(abc.ABC):
    """Base class for reaction-based queries."""

    @abc.abstractmethod
    def json(self) -> str:
        """Returns a JSON representation of the query."""

    @abc.abstractmethod
    def validate(self) -> None:
        """Checks the query for correctness.

        Raises:
            QueryException if the query is not valid.
        """

    @abc.abstractmethod
    def run(self, cursor: psycopg2.extensions.cursor, limit: Optional[int] = None) -> List[Result]:
        """Runs the query.

        Args:
            cursor: psycopg.cursor instance.
            limit: Integer maximum number of matches. If None (the default), no
                limit is set.

        Returns:
            List of Result instances.
        """


class RandomSampleQuery(ReactionQueryBase):
    """Takes a random sample of reactions."""

    def __init__(self, num_rows: int) -> None:
        """Initializes the query.

        Args:
            num_rows: Number of rows to return.
        """
        self._num_rows = num_rows

    def json(self) -> str:
        """Returns a JSON representation of the query."""
        return json.dumps({"num_rows": self._num_rows})

    def validate(self) -> None:
        """Checks the query for correctness.

        Raises:
            QueryException if the query is not valid.
        """
        if self._num_rows <= 0:
            raise QueryException("num_rows must be greater than zero")

    def run(self, cursor: psycopg2.extensions.cursor, limit: Optional[int] = None) -> List[Result]:
        """Runs the query.

        Args:
            cursor: psycopg.cursor instance.
            limit: Maximum number of matches. If None, no limit is set.

        Returns:
            List of Result instances.
        """
        del limit  # Unused.
        query = sql.SQL(
            """
            SELECT DISTINCT dataset.dataset_id, reaction.reaction_id, reaction.proto
            FROM ord.reaction TABLESAMPLE SYSTEM_ROWS (%s)
            JOIN dataset ON dataset.id = reaction.dataset_id
            """
        )
        args = [self._num_rows]
        logger.info("Running SQL command:%s", cursor.mogrify(query.as_string(cursor.connection), args).decode())
        cursor.execute(query, args)
        return fetch_results(cursor)


class DatasetIdQuery(ReactionQueryBase):
    """Looks up reactions by dataset ID."""

    def __init__(self, dataset_ids: List[str]) -> None:
        """Initializes the query.

        Args:
            dataset_ids: List of dataset IDs.
        """
        self._dataset_ids = dataset_ids

    def json(self) -> str:
        """Returns a JSON representation of the query."""
        return json.dumps({"datasetIds": self._dataset_ids})

    def validate(self) -> None:
        """Checks the query for correctness.

        Raises:
            QueryException if the query is not valid.
        """
        for dataset_id in self._dataset_ids:
            if not validations.is_valid_dataset_id(dataset_id):
                raise QueryException(f"invalid dataset ID: {dataset_id}")

    def run(self, cursor: psycopg2.extensions.cursor, limit: Optional[int] = None) -> List[Result]:
        """Runs the query.

        Args:
            cursor: psycopg.cursor instance.
            limit: Integer maximum number of matches. If None (the default), no
                limit is set.

        Returns:
            List of Result instances.
        """
        components = [
            sql.SQL(
                """
            SELECT DISTINCT dataset.dataset_id, reaction.reaction_id, reaction.proto
            FROM ord.reaction
            JOIN dataset ON dataset.id = reaction.dataset_id
            WHERE dataset.dataset_id = ANY (%s)"""
            )
        ]
        args = [self._dataset_ids]
        if limit:
            components.append(sql.SQL(" LIMIT %s"))
            args.append(limit)
        query = sql.Composed(components).join("")
        logger.info("Running SQL command:%s", cursor.mogrify(query.as_string(cursor.connection), args).decode())
        cursor.execute(query, args)
        return fetch_results(cursor)


class ReactionIdQuery(ReactionQueryBase):
    """Looks up reactions by ID."""

    def __init__(self, reaction_ids: List[str]) -> None:
        """Initializes the query.

        Args:
            reaction_ids: List of reaction IDs.
        """
        self._reaction_ids = reaction_ids

    def json(self) -> str:
        """Returns a JSON representation of the query."""
        return json.dumps({"reactionIds": self._reaction_ids})

    def validate(self) -> None:
        """Checks the query for correctness.

        Raises:
            QueryException if the query is not valid.
        """
        for reaction_id in self._reaction_ids:
            if not validations.is_valid_reaction_id(reaction_id):
                raise QueryException(f"invalid reaction ID: {reaction_id}")

    def run(self, cursor: psycopg2.extensions.cursor, limit: Optional[int] = None) -> List[Result]:
        """Runs the query.

        Args:
            cursor: psycopg.cursor instance.
            limit: Not used; present for compatibility.

        Returns:
            List of Result instances.
        """
        del limit  # Unused.
        query = sql.SQL(
            """
            SELECT DISTINCT dataset.dataset_id, reaction.reaction_id, reaction.proto
            FROM ord.reaction
            JOIN dataset ON dataset.id = reaction.dataset_id
            WHERE reaction.reaction_id = ANY (%s)
            """
        )
        args = [self._reaction_ids]
        logger.info("Running SQL command:%s", cursor.mogrify(query.as_string(cursor.connection), args).decode())
        cursor.execute(query, args)
        return fetch_results(cursor)


class ReactionSmartsQuery(ReactionQueryBase):
    """Matches reactions by reaction SMARTS."""

    def __init__(self, reaction_smarts: str) -> None:
        """Initializes the query.

        Args:
            reaction_smarts: Reaction SMARTS.
        """
        self._reaction_smarts = reaction_smarts

    def json(self) -> str:
        """Returns a JSON representation of the query."""
        return json.dumps({"reactionSmarts": self._reaction_smarts})

    def validate(self) -> None:
        """Checks the query for correctness.

        Raises:
            QueryException if the query is not valid.
        """
        try:
            rxn = rdChemReactions.ReactionFromSmarts(self._reaction_smarts)
            if not rxn:
                raise ValueError("cannot parse reaction SMARTS")
        except ValueError as error:
            raise QueryException(f"cannot parse reaction SMARTS: {self._reaction_smarts}") from error

    def run(self, cursor: psycopg2.extensions.cursor, limit: Optional[int] = None) -> List[Result]:
        """Runs the query.

        Args:
            cursor: psycopg.cursor instance.
            limit: Integer maximum number of matches. If None (the default), no
                limit is set.

        Returns:
            List of Result instances.
        """
        components = [
            sql.SQL(
                """
                SELECT DISTINCT dataset.dataset_id, reaction.reaction_id, reaction.proto
                FROM reaction
                JOIN rdkit.reactions ON rdkit.reactions.id = reaction.rdkit_reaction_id
                JOIN dataset ON dataset.id = reaction.dataset_id
                WHERE rdkit.reactions.reaction @> reaction_from_smarts(%s)
                """
            )
        ]
        args = [self._reaction_smarts]
        if limit:
            components.append(sql.SQL(" LIMIT %s"))
            args.append(limit)
        query = sql.Composed(components).join("")
        logger.info("Running SQL command:%s", cursor.mogrify(query.as_string(cursor.connection), args).decode())
        cursor.execute(query, args)
        return fetch_results(cursor)


class ReactionConversionQuery(ReactionQueryBase):
    """Looks up reactions by conversion."""

    def __init__(self, min_conversion: float, max_conversion: float) -> None:
        """Initializes the query.

        Args:
            min_conversion: Minimum conversion, as a percentage.
            max_conversion: Maximum conversion, as a percentage.
        """
        self._min_conversion = min_conversion
        self._max_conversion = max_conversion

    def json(self) -> str:
        """Returns a JSON representation of the query."""
        return json.dumps({"min_conversion": self._min_conversion, "max_conversion": self._max_conversion})

    def validate(self) -> None:
        """Checks the query for correctness.

        Raises:
            QueryException if the query is not valid.
        """
        # NOTE(skearnes): Reported values may be outside of [0, 100].

    def run(self, cursor: psycopg2.extensions.cursor, limit: Optional[int] = None) -> List[Result]:
        query = """
            SELECT DISTINCT dataset.dataset_id, reaction.reaction_id, reaction.proto
            FROM ord.reaction
            JOIN dataset ON dataset.id = reaction.dataset_id
            JOIN ord.reaction_outcome on reaction_outcome.reaction_id = reaction.id
            JOIN percentage on percentage.reaction_outcome_id = reaction_outcome.id
            WHERE percentage.value >= %s
              AND percentage.value <= %s
        """
        args = [self._min_conversion, self._max_conversion]
        if limit:
            query += "LIMIT %s"
            args.append(limit)
        logger.info("Running SQL command:%s", cursor.mogrify(query, args).decode())
        cursor.execute(query, args)
        return fetch_results(cursor)


class ReactionYieldQuery(ReactionQueryBase):
    """Looks up reactions by yield."""

    def __init__(self, min_yield: float, max_yield: float) -> None:
        """Initializes the query.

        Args:
            min_yield: Minimum yield, as a percentage.
            max_yield: Maximum yield, as a percentage.
        """
        self._min_yield = min_yield
        self._max_yield = max_yield

    def json(self) -> str:
        """Returns a JSON representation of the query."""
        return json.dumps({"min_yield": self._min_yield, "max_yield": self._max_yield})

    def validate(self) -> None:
        """Checks the query for correctness.

        Raises:
            QueryException if the query is not valid.
        """
        # NOTE(skearnes): Reported values may be outside of [0, 100].

    def run(self, cursor: psycopg2.extensions.cursor, limit: Optional[int] = None) -> List[Result]:
        query = """
            SELECT DISTINCT dataset.dataset_id, reaction.reaction_id, reaction.proto
            FROM ord.reaction
            JOIN dataset ON dataset.id = reaction.dataset_id
            JOIN ord.reaction_outcome on reaction_outcome.reaction_id = reaction.id
            JOIN ord.product_compound on product_compound.reaction_outcome_id = reaction_outcome.id
            JOIN ord.product_measurement on product_measurement.product_compound_id = product_compound.id
            JOIN percentage on percentage.product_measurement_id = product_measurement.id
            WHERE product_measurement.type = 'YIELD'
              AND percentage.value >= %s
              AND percentage.value <= %s
        """
        args = [self._min_yield, self._max_yield]
        if limit:
            query += "LIMIT %s"
            args.append(limit)
        logger.info("Running SQL command:%s", cursor.mogrify(query, args).decode())
        cursor.execute(query, args)
        return fetch_results(cursor)


class DoiQuery(ReactionQueryBase):
    """Looks up reactions by DOI."""

    def __init__(self, dois: List[str]) -> None:
        """Initializes the query.

        Args:
            dois: List of DOIs.
        """
        self._dois = dois

    def json(self) -> str:
        """Returns a JSON representation of the query."""
        return json.dumps({"dois": self._dois})

    def validate(self) -> None:
        """Checks the query for correctness.

        Raises:
            QueryException if the query is not valid.
        """
        for i, doi in enumerate(self._dois):
            try:
                parsed = message_helpers.parse_doi(doi)
            except ValueError as error:
                raise QueryException(f"invalid DOI: {doi}") from error
            if doi != parsed:
                # Trim the DOI as needed to match the database contents.
                logger.info(f"Updating DOI: {doi} -> {parsed}")
                self._dois[i] = parsed

    def run(self, cursor: psycopg2.extensions.cursor, limit: Optional[int] = None) -> List[Result]:
        """Runs the query.

        Args:
            cursor: psycopg2 cursor.
            limit: Integer maximum number of matches. If None (the default), no
                limit is set.

        Returns:
            List of Result instances.
        """
        components = [
            sql.SQL(
                """
                SELECT DISTINCT dataset.dataset_id, reaction.reaction_id, reaction.proto
                FROM ord.reaction
                JOIN dataset ON dataset.id = reaction.dataset_id
                JOIN reaction_provenance ON reaction_provenance.reaction_id = reaction.id
                WHERE reaction_provenance.doi = ANY (%s)
                """
            )
        ]
        args = [self._dois]
        if limit:
            components.append(sql.SQL(" LIMIT %s"))
            args.append(limit)
        query = sql.Composed(components).join("")
        logger.info("Running SQL command:%s", cursor.mogrify(query.as_string(cursor.connection), args).decode())
        cursor.execute(query, args)
        return fetch_results(cursor)


class ReactionComponentQuery(ReactionQueryBase):
    """Matches reactions by reaction component predicates."""

    def __init__(
        self, predicates: List[ReactionComponentPredicate], do_chiral_sss: bool = False, tanimoto_threshold: float = 0.5
    ) -> None:
        """Initializes the query.

        Args:
            predicates: List of ReactionComponentPredicate objects.
            do_chiral_sss: If True, consider stereochemistry in substructure searches.
            tanimoto_threshold: Float Tanimoto similarity threshold. Pairs below this threshold will not be
                considered matches.
        """
        self._predicates = predicates
        self._do_chiral_sss = do_chiral_sss
        self._tanimoto_threshold = tanimoto_threshold

    def json(self) -> str:
        """Returns a JSON representation of the query."""
        return json.dumps(
            {
                "useStereochemistry": self._do_chiral_sss,
                "similarity": self._tanimoto_threshold,
                "components": [predicate.to_dict() for predicate in self._predicates],
            }
        )

    def _setup(self, predicates: List[ReactionComponentPredicate], cursor: psycopg2.extensions.cursor) -> None:
        """Prepares the database for a query.

        Args:
            cursor: psycopg.cursor instance.
            predicates: Predicates included in this query.
        """
        command = sql.SQL("SET rdkit.do_chiral_sss=%s")
        args = [self._do_chiral_sss]
        logger.info("Running SQL command: %s", cursor.mogrify(command.as_string(cursor.connection), args).decode())
        cursor.execute(command, args)
        command = sql.SQL("SET rdkit.tanimoto_threshold=%s")
        tanimoto_threshold = self._tanimoto_threshold
        for predicate in predicates:
            if predicate.mode == ReactionComponentPredicate.MatchMode.EXACT:
                tanimoto_threshold = 1.0
        args = [tanimoto_threshold]
        logger.info("Running SQL command: %s", cursor.mogrify(command.as_string(cursor.connection), args).decode())
        cursor.execute(command, args)

    def validate(self) -> None:
        """Checks the query for correctness.

        Raises:
            QueryException if the query is not valid.
        """
        for predicate in self._predicates:
            if predicate.mode == predicate.MatchMode.SMARTS:
                mol = Chem.MolFromSmarts(predicate.pattern)
            else:
                mol = Chem.MolFromSmiles(predicate.pattern)
            if not mol:
                raise QueryException(f"cannot parse pattern: {predicate.pattern}")

    def run(self, cursor: psycopg2.extensions.cursor, limit: Optional[int] = None) -> List[Result]:
        """Runs the query.

        Args:
            cursor: psycopg.cursor instance.
            limit: Integer maximum number of matches. If None (the default), no
                limit is set.

        Returns:
            List of Result instances.
        """
        # The RDKit Postgres cartridge only allows the Tanimoto threshold to be
        # set once per query. If we mix EXACT and SIMILAR modes, we need to get
        # the results in separate queries and compute the intersection.
        exact, other = [], []
        for predicate in self._predicates:
            if predicate.mode == ReactionComponentPredicate.MatchMode.EXACT:
                exact.append(predicate)
            else:
                other.append(predicate)
        exact_results = self._run(exact, cursor=cursor, limit=limit)
        other_results = self._run(other, cursor=cursor, limit=limit)
        if exact_results and other_results:
            return list(set(exact_results).intersection(set(other_results)))
        if exact_results:
            return exact_results
        if other_results:
            return other_results
        return []

    def _run(
        self,
        predicates: List[ReactionComponentPredicate],
        cursor: psycopg2.extensions.cursor,
        limit: Optional[int] = None,
    ) -> List[Result]:
        """Runs the query for a set of predicates."""
        if not predicates:
            return []
        self._setup(predicates, cursor)
        components = []
        args = []
        for predicate in predicates:
            if predicate.target == ReactionComponentPredicate.Target.INPUT:
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
            predicate_sql, predicate_args = predicate.get()
            components.append(
                f"""
                SELECT DISTINCT dataset.dataset_id, reaction.reaction_id, reaction.proto
                FROM ord.reaction
                {mols_sql}
                JOIN ord.dataset ON dataset.id = reaction.dataset_id
                WHERE
                {predicate_sql}
                """
            )
            args.extend(predicate_args)
        assert len(components) > 0
        query = "\nINTERSECT\n".join(components)
        if limit:
            query += "\nLIMIT %s\n"
            args.append(limit)
        try:
            logger.info("Running SQL command:%s", cursor.mogrify(query, args).decode())
        except IndexError as error:
            raise ValueError((query, args)) from error
        cursor.execute(query, args)
        return fetch_results(cursor)


class ReactionComponentPredicate:
    """Specifies a single reaction component predicate."""

    class Target(enum.Enum):
        """Search targets."""

        INPUT = enum.auto()
        OUTPUT = enum.auto()

        @classmethod
        def from_name(cls, name: str) -> ReactionComponentPredicate.Target:
            """Takes a matching criterion from a URL param."""
            return cls[name.upper()]

    class MatchMode(enum.Enum):
        """Interpretations for SMILES and SMARTS strings."""

        EXACT = enum.auto()
        SIMILAR = enum.auto()
        SUBSTRUCTURE = enum.auto()
        SMARTS = enum.auto()

        @classmethod
        def from_name(cls, name: str) -> ReactionComponentPredicate.MatchMode:
            """Takes a matching criterion from a URL param."""
            return cls[name.upper()]

    def __init__(self, pattern: str, target: Target, mode: MatchMode) -> None:
        """Initializes the predicate.

        Args:
            pattern: SMILES or SMARTS pattern.
            target: Search target.
            mode: ReactionComponentPredicate.MatchMode.

        Raises:
            ValueError: If `table` is not allowed.
        """
        self._pattern = pattern
        self._target = target
        self._mode = mode

    @property
    def pattern(self) -> str:
        return self._pattern

    @property
    def target(self) -> Target:
        return self._target

    @property
    def mode(self) -> MatchMode:
        return self._mode

    def to_dict(self) -> Dict[str, str]:
        """Returns a dict representation of the predicate."""
        return {"pattern": self._pattern, "target": self._target.name.lower(), "mode": self._mode.name.lower()}

    def get(self) -> Tuple[str, List[str]]:
        """Builds the SQL predicate.

        Returns:
            predicate: sql.SQL query object.
            args: List of arguments for `predicate`.
        """
        if self._mode in [ReactionComponentPredicate.MatchMode.SIMILAR, ReactionComponentPredicate.MatchMode.EXACT]:
            predicate = "rdkit.mols.morgan_bfp %% morganbv_fp(%s)"  # Escape the % operator.
        elif self._mode == ReactionComponentPredicate.MatchMode.SUBSTRUCTURE:
            predicate = "rdkit.mols.mol @> %s"
        elif self._mode == ReactionComponentPredicate.MatchMode.SMARTS:
            predicate = "rdkit.mols.mol @> %s::qmol"
        else:
            raise ValueError(f"unsupported mode: {self._mode}")
        return predicate, [self._pattern]


class OrdPostgres:
    """Class for performing SQL queries on the ORD."""

    def __init__(self, dbname: str, user: str, password: str, host: str, port: int) -> None:
        """Initializes an instance of OrdPostgres.

        Args:
            dbname: Text database name.
            user: Text user name.
            password: Text user password.
            host: Text host name.
            port: Integer port.
        """
        self._connection = psycopg2.connect(
            dbname=dbname, user=user, password=password, host=host, port=port, options="-c search_path=public,ord"
        )
        self._connection.set_session(readonly=True)

    @property
    def connection(self) -> psycopg2.extensions.connection:
        return self._connection

    def cursor(self) -> psycopg2.extensions.cursor:
        return self._connection.cursor()

    def run_query(
        self, query: ReactionQueryBase, limit: Optional[int] = None, return_ids: bool = False
    ) -> List[Result]:
        """Runs a query against the database.

        Args:
            query: ReactionQueryBase query.
            limit: Integer maximum number of matches. If None (the default), no
                limit is set.
            return_ids: If True, only return reaction IDs. If False, return
                full Reaction records.

        Returns:
            List of Result instances.
        """
        query.validate()
        with self._connection, self.cursor() as cursor:
            results = query.run(cursor, limit=limit)
            self._connection.rollback()  # Revert rdkit runtime configuration.
        if return_ids:
            only_ids = []
            for result in results:
                only_ids.append(Result(dataset_id=result.dataset_id, reaction_id=result.reaction_id))
            return only_ids
        return results


class QueryException(Exception):
    """Exception class for ORD queries."""
