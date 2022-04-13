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
"""Tests for ord_interface.query."""

from absl.testing import absltest
from absl.testing import parameterized
import numpy as np

import ord_interface
from ord_interface import query


class QueryTest(parameterized.TestCase, absltest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.postgres = query.OrdPostgres(
            dbname=ord_interface.POSTGRES_DB,
            user=ord_interface.POSTGRES_USER,
            password=ord_interface.POSTGRES_PASSWORD,
            # Matches the service name in docker-compose.yml.
            host='ord-postgres',
            port=ord_interface.POSTGRES_PORT)

    def test_random_sample_query(self):
        command = query.RandomSampleQuery(16)
        results = self.postgres.run_query(command, return_ids=True)
        self.assertLen(results, 16)

    def test_dataset_id_query(self):
        dataset_ids = ['ord_dataset-46ff9a32d9e04016b9380b1b1ef949c3']
        command = query.DatasetIdQuery(dataset_ids)
        results = self.postgres.run_query(command, limit=10, return_ids=True)
        self.assertLen(results, 10)

    def test_reaction_id_query(self):
        reaction_ids = ['ord-e49ed67da61e4cddabd3c84a72fed227']
        command = query.ReactionIdQuery(reaction_ids)
        results = self.postgres.run_query(command, limit=10, return_ids=False)
        self.assertCountEqual([result.reaction_id for result in results],
                              reaction_ids)

    def test_reaction_smarts_query(self):
        pattern = '[#6]>>[#7]'
        command = query.ReactionSmartsQuery(pattern)
        results = self.postgres.run_query(command, limit=10, return_ids=True)
        self.assertLen(results, 10)

    def test_doi_query(self):
        dois = ['10.1021/acscatal.0c02247']
        command = query.DoiQuery(dois)
        results = self.postgres.run_query(command, limit=10)
        self.assertLen(results, 10)
        for result in results:
            self.assertIn(result.reaction.provenance.doi, dois)

    def test_exact_query(self):
        pattern = 'O=C1NC2=CC(Br)=CC3=C2N(C(CC(OC)=O)CC3)C1=O'
        mode = query.ReactionComponentPredicate.MatchMode.EXACT
        predicates = [
            query.ReactionComponentPredicate(pattern, table='inputs', mode=mode)
        ]
        command = query.ReactionComponentQuery(predicates)
        results = self.postgres.run_query(command, limit=5, return_ids=True)
        self.assertLen(results, 5)
        # Check that we remove redundant reaction IDs.
        reaction_ids = [result.reaction_id for result in results]
        self.assertCountEqual(reaction_ids, np.unique(reaction_ids))

    def test_substructure_query(self):
        pattern = 'C'
        mode = query.ReactionComponentPredicate.MatchMode.SUBSTRUCTURE
        predicates = [
            query.ReactionComponentPredicate(pattern, table='inputs', mode=mode)
        ]
        command = query.ReactionComponentQuery(predicates)
        results = self.postgres.run_query(command, limit=10, return_ids=True)
        self.assertLen(results, 10)
        # Check that we remove redundant reaction IDs.
        reaction_ids = [result.reaction_id for result in results]
        self.assertCountEqual(reaction_ids, np.unique(reaction_ids))

    def test_smarts_query(self):
        pattern = '[#6]'
        mode = query.ReactionComponentPredicate.MatchMode.SMARTS
        predicates = [
            query.ReactionComponentPredicate(pattern, table='inputs', mode=mode)
        ]
        command = query.ReactionComponentQuery(predicates)
        results = self.postgres.run_query(command, limit=10, return_ids=True)
        self.assertLen(results, 10)
        # Check that we remove redundant reaction IDs.
        reaction_ids = [result.reaction_id for result in results]
        self.assertCountEqual(reaction_ids, np.unique(reaction_ids))

    def test_similarity_query(self):
        pattern = 'CC=O'
        mode = query.ReactionComponentPredicate.MatchMode.SIMILAR
        predicates = [
            query.ReactionComponentPredicate(pattern, table='inputs', mode=mode)
        ]
        command = query.ReactionComponentQuery(predicates,
                                               tanimoto_threshold=0.5)
        results = self.postgres.run_query(command, limit=10, return_ids=True)
        self.assertEmpty(results)
        command = query.ReactionComponentQuery(predicates,
                                               tanimoto_threshold=0.05)
        results = self.postgres.run_query(command, limit=10, return_ids=True)
        self.assertLen(results, 10)
        # Check that we remove redundant reaction IDs.
        reaction_ids = [result.reaction_id for result in results]
        self.assertCountEqual(reaction_ids, np.unique(reaction_ids))

    def test_bad_smiles(self):
        pattern = 'invalid_smiles'
        mode = query.ReactionComponentPredicate.MatchMode.SUBSTRUCTURE
        predicates = [
            query.ReactionComponentPredicate(pattern, table='inputs', mode=mode)
        ]
        command = query.ReactionComponentQuery(predicates)
        with self.assertRaisesRegex(query.QueryException,
                                    'cannot parse pattern: invalid_smiles'):
            self.postgres.run_query(command)


if __name__ == '__main__':
    absltest.main()
