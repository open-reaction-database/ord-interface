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
"""Tests for ord_interface.ord_client."""

from absl.testing import absltest
from absl.testing import parameterized

from ord_interface.client import ord_client


class OrdClientTest(parameterized.TestCase, absltest.TestCase):

    def setUp(self):
        super().setUp()
        self.client = ord_client.OrdClient(target='http://localhost:5000')

    @parameterized.parameters(
        ('ord_dataset-d319c2a22ecf4ce59db1a18ae71d529c', 264))
    def test_fetch_dataset(self, dataset_id, expected_num_reactions):
        dataset = self.client.fetch_dataset(dataset_id)
        self.assertLen(dataset.reactions, expected_num_reactions)

    @parameterized.parameters(
        (['ord_dataset-d319c2a22ecf4ce59db1a18ae71d529c'], [264]))
    def test_fetch_datasets(self, dataset_ids, expected_num_reactions):
        datasets = self.client.fetch_datasets(dataset_ids)
        self.assertLen(dataset_ids, len(expected_num_reactions))
        self.assertLen(datasets, len(expected_num_reactions))
        for dataset, expected in zip(datasets, expected_num_reactions):
            self.assertLen(dataset.reactions, expected)

    @parameterized.parameters(
        ('ord-cf0d04017ede4c8aab8a15119c53e57b', 'Connor W. Coley'),)
    def test_fetch_reaction(self, reaction_id, created_by):
        reaction = self.client.fetch_reaction(reaction_id)
        self.assertEqual(reaction.provenance.record_created.person.name,
                         created_by)

    @parameterized.parameters(
        ([
            'ord-cf0d04017ede4c8aab8a15119c53e57b',
            'ord-153de2a96e07484e93d525d81f966789'
        ], ['Connor W. Coley', 'Connor W. Coley']),)
    def test_fetch_reactions(self, reaction_ids, created_by):
        reactions = self.client.fetch_reactions(reaction_ids)
        self.assertLen(reaction_ids, len(created_by))
        self.assertLen(reactions, len(created_by))
        for reaction, expected in zip(reactions, created_by):
            self.assertEqual(reaction.provenance.record_created.person.name,
                             expected)

    def test_query_dataset_ids(self):
        results = self.client.query(
            dataset_ids=['ord_dataset-89b083710e2d441aa0040c361d63359f'])
        self.assertLen(results, 24)

    def test_query_reaction_ids(self):
        results = self.client.query(reaction_ids=[
            'ord-cf0d04017ede4c8aab8a15119c53e57b',
            'ord-153de2a96e07484e93d525d81f966789',
        ])
        self.assertLen(results, 2)

    def test_query_reaction_smarts(self):
        results = self.client.query(
            reaction_smarts='[Br]C1=CC=C(C(C)=O)C=C1>CN(C)C=O>')
        self.assertLen(results, 9)

    def test_query_dois(self):
        doi = '10.1126/science.1255525'
        results = self.client.query(dois=[doi])
        self.assertLen(results, 24)
        for result in results:
            self.assertEqual(result.reaction.provenance.doi, doi)

    def test_query_single_component(self):
        component = ord_client.ComponentQuery('[Br]C1=CC=C(C(C)=O)C=C1',
                                              source='input',
                                              mode='exact')
        results = self.client.query(components=[component])
        self.assertLen(results, 10)

    def test_query_multiple_components(self):
        component1 = ord_client.ComponentQuery('[Br]C1=CC=C(C(C)=O)C=C1',
                                               source='input',
                                               mode='exact')
        component2 = ord_client.ComponentQuery('CC(C)(C)OC(=O)N1CCCCC1C(=O)O',
                                               source='input',
                                               mode='exact')
        results = self.client.query(components=[component1, component2])
        self.assertLen(results, 1)

    def test_query_stereochemistry(self):
        # TODO(kearnes): Add a test once we have more chiral stuff.
        pass

    def test_query_similarity(self):
        component = ord_client.ComponentQuery('[Br]C1=CC=C(C(C)=O)C=C1',
                                              source='input',
                                              mode='similar')
        results = self.client.query(components=[component], similarity=0.6)
        self.assertLen(results, 10)


if __name__ == '__main__':
    absltest.main()
