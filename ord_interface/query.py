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

import abc
import binascii
import enum
import json
from typing import Dict, List, Optional

from absl import logging
import psycopg2
import psycopg2.extensions
from psycopg2 import sql
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2


def fetch_results(cursor):
    """Fetches query results.

    Args:
        cursor: psycopg.cursor instance.

    Returns:
        Dict mapping reaction IDs to serialized Reaction protos.
    """
    reactions = {}
    for reaction_id, serialized in cursor:
        reactions[reaction_id] = serialized
    return reactions


class ReactionQueryBase(abc.ABC):
    """Base class for reaction-based queries."""

    @abc.abstractmethod
    def json(self):
        """Returns a JSON representation of the query."""

    @abc.abstractmethod
    def validate(self):
        """Checks the query for correctness.

        Raises:
            QueryException if the query is not valid.
        """

    @abc.abstractmethod
    def run(self, cursor, limit=None):
        """Runs the query.

        Args:
            cursor: psycopg.cursor instance.
            limit: Integer maximum number of matches. If None (the default), no
                limit is set.

        Returns:
            Dict mapping reaction IDs to serialized Reaction protos.
        """


class RandomSampleQuery(ReactionQueryBase):
    """Takes a random sample of reactions."""

    def __init__(self, num_rows: int):
        """Initializes the query.

        Args:
            num_rows: Number of rows to return.
        """
        self._num_rows = num_rows

    def json(self):
        """Returns a JSON representation of the query."""
        return json.dumps({'num_rows': self._num_rows})

    def validate(self):
        """Checks the query for correctness.

        Raises:
            QueryException if the query is not valid.
        """
        if self._num_rows <= 0:
            raise QueryException('num_rows must be greater than zero')

    def run(self,
            cursor: psycopg2.extensions.cursor,
            limit: Optional[int] = None):
        """Runs the query.

        Args:
            cursor: psycopg.cursor instance.
            limit: Maximum number of matches. If None, no limit is set.

        Returns:
            Dict mapping reaction IDs to serialized Reaction protos.
        """
        del limit  # Unused.
        query = sql.SQL("""
            SELECT DISTINCT reaction_id, serialized 
            FROM reactions 
            TABLESAMPLE SYSTEM_ROWS (%s)""")
        args = [self._num_rows]
        logging.info('Running SQL command:%s',
                     cursor.mogrify(query, args).decode())
        cursor.execute(query, args)
        return fetch_results(cursor)


class DatasetIdQuery(ReactionQueryBase):
    """Looks up reactions by dataset ID."""

    def __init__(self, dataset_ids):
        """Initializes the query.

        Args:
            dataset_ids: List of dataset IDs.
        """
        self._dataset_ids = dataset_ids

    def json(self):
        """Returns a JSON representation of the query."""
        return json.dumps({'datasetIds': self._dataset_ids})

    def validate(self):
        """Checks the query for correctness.

        Raises:
            QueryException if the query is not valid.
        """
        for dataset_id in self._dataset_ids:
            if not validations.is_valid_dataset_id(dataset_id):
                raise QueryException(f'invalid dataset ID: {dataset_id}')

    def run(self, cursor, limit=None):
        """Runs the query.

        Args:
            cursor: psycopg.cursor instance.
            limit: Integer maximum number of matches. If None (the default), no
                limit is set.

        Returns:
            Dict mapping reaction IDs to serialized Reaction protos.
        """
        components = [
            sql.SQL("""
            SELECT DISTINCT reaction_id, serialized 
            FROM reactions 
            WHERE dataset_id = ANY (%s)""")
        ]
        args = [self._dataset_ids]
        if limit:
            components.append(sql.SQL(' LIMIT %s'))
            args.append(limit)
        query = sql.Composed(components).join('')
        logging.info('Running SQL command:%s',
                     cursor.mogrify(query, args).decode())
        cursor.execute(query, args)
        return fetch_results(cursor)


class ReactionIdQuery(ReactionQueryBase):
    """Looks up reactions by ID."""

    def __init__(self, reaction_ids):
        """Initializes the query.

        Args:
            reaction_ids: List of reaction IDs.
        """
        self._reaction_ids = reaction_ids

    def json(self):
        """Returns a JSON representation of the query."""
        return json.dumps({'reactionIds': self._reaction_ids})

    def validate(self):
        """Checks the query for correctness.

        Raises:
            QueryException if the query is not valid.
        """
        for reaction_id in self._reaction_ids:
            if not validations.is_valid_reaction_id(reaction_id):
                raise QueryException(f'invalid reaction ID: {reaction_id}')

    def run(self, cursor, limit=None):
        """Runs the query.

        Args:
            cursor: psycopg.cursor instance.
            limit: Not used; present for compatibility.

        Returns:
            Dict mapping reaction IDs to serialized Reaction protos.
        """
        del limit  # Unused.
        query = sql.SQL("""
            SELECT DISTINCT reaction_id, serialized 
            FROM reactions 
            WHERE reaction_id = ANY (%s)""")
        args = [self._reaction_ids]
        logging.info('Running SQL command:%s',
                     cursor.mogrify(query, args).decode())
        cursor.execute(query, args)
        return fetch_results(cursor)


class ReactionSmartsQuery(ReactionQueryBase):
    """Matches reactions by reaction SMARTS."""

    def __init__(self, reaction_smarts):
        """Initializes the query.

        Args:
            reaction_smarts: Reaction SMARTS.
        """
        self._reaction_smarts = reaction_smarts

    def json(self):
        """Returns a JSON representation of the query."""
        return json.dumps({'reactionSmarts': self._reaction_smarts})

    def validate(self):
        """Checks the query for correctness.

        Raises:
            QueryException if the query is not valid.
        """
        try:
            rxn = rdChemReactions.ReactionFromSmarts(self._reaction_smarts)
            if not rxn:
                raise ValueError('cannot parse reaction SMARTS')
        except ValueError as error:
            raise QueryException(
                f'cannot parse reaction SMARTS: {self._reaction_smarts}'
            ) from error

    def run(self, cursor, limit=None):
        """Runs the query.

        Args:
            cursor: psycopg.cursor instance.
            limit: Integer maximum number of matches. If None (the default), no
                limit is set.

        Returns:
            Dict mapping reaction IDs to serialized Reaction protos.
        """
        components = [
            sql.SQL("""
            SELECT DISTINCT reaction_id, serialized
            FROM reactions
            INNER JOIN rdk.reactions USING (reaction_id)
            WHERE rdk.reactions.r@>reaction_from_smarts(%s)""")
        ]
        args = [self._reaction_smarts]
        if limit:
            components.append(sql.SQL(' LIMIT %s'))
            args.append(limit)
        query = sql.Composed(components).join('')
        logging.info('Running SQL command:%s',
                     cursor.mogrify(query, args).decode())
        cursor.execute(query, args)
        return fetch_results(cursor)


class DoiQuery(ReactionQueryBase):
    """Looks up reactions by DOI."""

    def __init__(self, dois: List[str]):
        """Initializes the query.

        Args:
            dois: List of DOIs.
        """
        self._dois = dois

    def json(self) -> str:
        """Returns a JSON representation of the query."""
        return json.dumps({'dois': self._dois})

    def validate(self):
        """Checks the query for correctness.

        Raises:
            QueryException if the query is not valid.
        """
        for i, doi in enumerate(self._dois):
            try:
                parsed = message_helpers.parse_doi(doi)
            except ValueError as error:
                raise QueryException(f'invalid DOI: {doi}') from error
            if doi != parsed:
                # Trim the DOI as needed to match the database contents.
                logging.info(f'Updating DOI: {doi} -> {parsed}')
                self._dois[i] = parsed

    def run(self,
            cursor: psycopg2.extensions.cursor,
            limit: Optional[int] = None) -> Dict[str, reaction_pb2.Reaction]:
        """Runs the query.

        Args:
            cursor: psycopg2 cursor.
            limit: Integer maximum number of matches. If None (the default), no
                limit is set.

        Returns:
            Dict mapping reaction IDs to serialized Reaction protos.
        """
        components = [
            sql.SQL("""
            SELECT DISTINCT reaction_id, serialized 
            FROM reactions 
            WHERE doi = ANY (%s)""")
        ]
        args = [self._dois]
        if limit:
            components.append(sql.SQL(' LIMIT %s'))
            args.append(limit)
        query = sql.Composed(components).join('')
        logging.info('Running SQL command:%s',
                     cursor.mogrify(query, args).decode())
        cursor.execute(query, args)
        return fetch_results(cursor)


class ReactionComponentQuery(ReactionQueryBase):
    """Matches reactions by reaction component predicates."""

    def __init__(self, predicates, do_chiral_sss=False, tanimoto_threshold=0.5):
        """Initializes the query.

        Args:
            predicates: List of ReactionComponentPredicate objects.
            do_chiral_sss: If True, consider stereochemistry in substructure
                searches.
            tanimoto_threshold: Float Tanimoto similarity threshold. Pairs
                below this threshold will not be considered matches.
        """
        self._predicates = predicates
        self._do_chiral_sss = do_chiral_sss
        self._tanimoto_threshold = tanimoto_threshold

    def json(self):
        """Returns a JSON representation of the query."""
        return json.dumps({
            'useStereochemistry':
                self._do_chiral_sss,
            'similarity':
                self._tanimoto_threshold,
            'components': [
                predicate.to_dict() for predicate in self._predicates
            ],
        })

    def _setup(self, cursor):
        """Prepares the database for a query.

        Args:
            cursor: psycopg.cursor instance.
        """
        command = sql.SQL('SET rdkit.do_chiral_sss=%s')
        args = [self._do_chiral_sss]
        logging.info('Running SQL command: %s',
                     cursor.mogrify(command, args).decode())
        cursor.execute(command, args)
        command = sql.SQL('SET rdkit.tanimoto_threshold=%s')
        args = [self._tanimoto_threshold]
        logging.info('Running SQL command: %s',
                     cursor.mogrify(command, args).decode())
        cursor.execute(command, args)

    def _get_tables(self):
        """Identifies the minimum set of tables to join for the query."""
        tables = []
        requires_inputs = False
        requires_rdk_inputs = False
        requires_outputs = False
        requires_rdk_outputs = False
        for predicate in self._predicates:
            if predicate.table == 'inputs':
                if predicate.mode == ReactionComponentPredicate.MatchMode.EXACT:
                    requires_inputs = True
                else:
                    requires_rdk_inputs = True
            elif predicate.table == 'outputs':
                if predicate.mode == ReactionComponentPredicate.MatchMode.EXACT:
                    requires_outputs = True
                else:
                    requires_rdk_outputs = True
        if requires_inputs:
            tables.append(
                sql.SQL("""
            INNER JOIN inputs USING (reaction_id) """))
        if requires_rdk_inputs:
            tables.append(
                sql.SQL("""
            INNER JOIN rdk.inputs USING (reaction_id) """))
        if requires_outputs:
            tables.append(
                sql.SQL("""
            INNER JOIN outputs USING (reaction_id) """))
        if requires_rdk_outputs:
            tables.append(
                sql.SQL("""
            INNER JOIN rdk.outputs USING (reaction_id) """))
        return tables

    def validate(self):
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
                raise QueryException(
                    f'cannot parse pattern: {predicate.pattern}')

    def run(self, cursor, limit=None):
        """Runs the query.

        Args:
            cursor: psycopg.cursor instance.
            limit: Integer maximum number of matches. If None (the default), no
                limit is set.

        Returns:
            Dict mapping reaction IDs to serialized Reaction protos.
        """
        if not self._predicates:
            return {}
        self._setup(cursor)
        predicates = []
        args = []
        for predicate in self._predicates:
            components = [
                sql.SQL("""
                SELECT DISTINCT reaction_id, serialized
                FROM reactions """)
            ]
            components.extend(self._get_tables())
            components.append(sql.SQL("""
                WHERE """))
            predicate_sql, predicate_args = predicate.get()
            components.append(predicate_sql)
            args.extend(predicate_args)
            predicates.append(sql.Composed(components))
        components = [sql.Composed(predicates).join(' INTERSECT ')]
        if limit:
            components.append(sql.SQL(' LIMIT %s'))
            args.append(limit)
        query = sql.Composed(components).join('')
        logging.info('Running SQL command:%s',
                     cursor.mogrify(query, args).decode())
        cursor.execute(query, args)
        return fetch_results(cursor)


class ReactionComponentPredicate:
    """Specifies a single reaction component predicate."""

    _ALLOWED_TABLES = ['inputs', 'outputs']
    SOURCE_TO_TABLE = {'input': 'inputs', 'output': 'outputs'}
    _TABLE_TO_SOURCE = {'inputs': 'input', 'outputs': 'output'}

    class MatchMode(enum.Enum):
        """Interpretations for SMILES and SMARTS strings."""
        EXACT = 1
        SIMILAR = 2
        SUBSTRUCTURE = 3
        SMARTS = 4

        @classmethod
        def from_name(cls, name):
            """Takes a matching criterion from a URL param."""
            return ReactionComponentPredicate.MatchMode[name.upper()]

    def __init__(self, pattern, table, mode):
        """Initializes the predicate.

        Args:
            pattern: SMILES or SMARTS pattern.
            table: Table to search.
            mode: ReactionComponentPredicate.MatchMode.

        Raises:
            ValueError: If `table` is not allowed.
        """
        if table not in self._ALLOWED_TABLES:
            raise ValueError(f'table must be in {self._ALLOWED_TABLES}')
        self._pattern = pattern
        self._table = table
        self._mode = mode

    @property
    def pattern(self):
        return self._pattern

    @property
    def table(self):
        return self._table

    @property
    def mode(self):
        return self._mode

    def to_dict(self):
        """Returns a dict representation of the predicate."""
        return {
            'pattern': self._pattern,
            'source': self._TABLE_TO_SOURCE[self._table],
            'mode': self._mode.name.lower(),
        }

    def get(self):
        """Builds the SQL predicate.

        Returns:
            predicate: sql.SQL query object.
            args: List of arguments for `predicate`.
        """
        if self._mode == ReactionComponentPredicate.MatchMode.EXACT:
            # Canonicalize the SMILES.
            self._pattern = Chem.MolToSmiles(Chem.MolFromSmiles(self._pattern))
            predicate = sql.SQL('{} = %s').format(
                sql.Identifier(self._table, 'smiles'))
        elif self._mode == ReactionComponentPredicate.MatchMode.SIMILAR:
            predicate = sql.SQL('{}%%morganbv_fp(%s)').format(
                sql.Identifier('rdk', self._table, 'mfp2'))
        elif self._mode == ReactionComponentPredicate.MatchMode.SUBSTRUCTURE:
            predicate = sql.SQL('{}@>%s').format(
                sql.Identifier('rdk', self._table, 'm'))
        elif self._mode == ReactionComponentPredicate.MatchMode.SMARTS:
            predicate = sql.SQL('{}@>%s::qmol').format(
                sql.Identifier('rdk', self._table, 'm'))
        else:
            # NOTE(kearnes): This should never happen, so I'm leaving it out of
            # the docstring.
            raise ValueError(f'unsupported mode: {self._mode}')
        return predicate, [self._pattern]


class OrdPostgres:
    """Class for performing SQL queries on the ORD."""

    def __init__(self, dbname, user, password, host, port):
        """Initializes an instance of OrdPostgres.

        Args:
            dbname: Text database name.
            user: Text user name.
            password: Text user password.
            host: Text host name.
            port: Integer port.
        """
        self._connection = psycopg2.connect(dbname=dbname,
                                            user=user,
                                            password=password,
                                            host=host,
                                            port=port)
        self._connection.set_session(readonly=True)

    def cursor(self):
        return self._connection.cursor()

    def run_query(self, query, limit=None, return_ids=False):
        """Runs a query against the database.

        Args:
            query: ReactionQueryBase query.
            limit: Integer maximum number of matches. If None (the default), no
                limit is set.
            return_ids: If True, only return reaction IDs. If False, return
                full Reaction records.

        Returns:
            dataset_pb2.Dataset containing the matched reactions (or IDs).
        """
        query.validate()
        with self._connection, self.cursor() as cursor:
            reactions = query.run(cursor, limit=limit)
            self._connection.rollback()  # Revert rdkit runtime configuration.
        if return_ids:
            return dataset_pb2.Dataset(reaction_ids=reactions.keys())
        unserialized = []
        for serialized in reactions.values():
            reaction = reaction_pb2.Reaction.FromString(
                binascii.unhexlify(serialized.tobytes()))
            unserialized.append(reaction)
        return dataset_pb2.Dataset(reactions=unserialized)


class QueryException(Exception):
    """Exception class for ORD queries."""
