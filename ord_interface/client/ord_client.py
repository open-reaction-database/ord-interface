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
"""Python API for the Open Reaction Database."""
from __future__ import annotations

import binascii
from typing import List, Optional

import requests

from ord_schema import message_helpers
from ord_schema import validations
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

from ord_interface.client import query

TARGET = "https://open-reaction-database.org"


def parse_response(response: requests.Response) -> list[query.Result]:
    """Parses a JSON response into Result objects."""
    results = []
    for result in response.json():
        result["proto"] = binascii.unhexlify(result["proto"])
        results.append(query.Result(**result))
    return results


class OrdClient:
    """Client for the Open Reaction Database."""

    def __init__(self, target: Optional[str] = None, prefix: str = "") -> None:
        """Initializes the client.

        Args:
            target: Endpoint URL for client queries. Defaults to the public ORD client URL (defined in TARGET).
            prefix: URL prefix.
        """
        if not target:
            target = TARGET
        self._target = target
        self._prefix = prefix

    def fetch_datasets(self, dataset_ids: List[str]) -> List[dataset_pb2.Dataset]:
        """Fetches one or more Dataset messages.

        Args:
            dataset_ids: List of dataset IDs.

        Returns:
            List of Dataset messages.
        """
        # NOTE(kearnes): Return datasets in the same order as requested.
        datasets = {dataset_id: self.fetch_dataset(dataset_id) for dataset_id in set(dataset_ids)}
        return [datasets[dataset_id] for dataset_id in dataset_ids]

    @staticmethod
    def fetch_dataset(dataset_id: str) -> dataset_pb2.Dataset:
        """Fetches a single Dataset message.

        Args:
            dataset_id: String dataset ID.

        Returns:
            Dataset message.
        """
        return message_helpers.fetch_dataset(dataset_id)

    def fetch_reactions(self, reaction_ids: List[str], timeout: float = 10.0) -> List[reaction_pb2.Reaction]:
        """Fetches one or more Reaction messages.

        Args:
            reaction_ids: List of reaction IDs.
            timeout: Number of seconds to wait before timing out the request.

        Returns:
            List of Reaction messages.

        Raises:
            ValueError: A reaction ID is invalid.
        """
        for reaction_id in reaction_ids:
            if not validations.is_valid_reaction_id(reaction_id):
                raise ValueError(f"Invalid reaction ID: {reaction_id}")
        target = self._target + self._prefix + "/api/fetch_reactions"
        response = requests.post(target, json=list(set(reaction_ids)), timeout=timeout)
        results = parse_response(response)
        reactions = {result.reaction_id: result.reaction for result in results}
        # Return reactions in the same order as requested.
        return [reactions[reaction_id] for reaction_id in reaction_ids]

    def fetch_reaction(self, reaction_id: str) -> reaction_pb2.Reaction:
        """Fetches a single Reaction message.

        Args:
            reaction_id: String reaction ID.

        Returns:
            Reaction message.
        """
        reactions = self.fetch_reactions([reaction_id])
        assert len(reactions) == 1
        return reactions[0]

    def query(  # pylint: disable=too-many-arguments,too-many-locals
        self,
        dataset_ids: Optional[List[str]] = None,
        reaction_ids: Optional[List[str]] = None,
        reaction_smarts: Optional[str] = None,
        dois: Optional[List[str]] = None,
        components: Optional[List[ComponentQuery]] = None,
        use_stereochemistry: Optional[bool] = None,
        similarity: Optional[float] = None,
        timeout: float = 10.0,
    ) -> List[query.Result]:
        """Executes a query against the Open Reaction Database.

        Args:
            dataset_ids: List of dataset IDs to fetch. This is provided for
                completeness, but fetch_dataset(s) should be preferred.
            reaction_ids: List of reaction IDs to fetch. This is provided for
                completeness, but fetch_reaction(s) should be preferred.
            reaction_smarts: Reaction SMARTS pattern.
            dois: List of DOIs to match.
            components: List of ComponentQuery instances.
            use_stereochemistry: Boolean whether to use stereochemistry when
                matching.
            similarity: Float similarity threshold for SIMILAR queries.
            timeout: Number of seconds to wait before timing out the request.

        Returns:
            List of Result instances.

        Raises:
            ValueError: A reaction ID is invalid.
        """
        if dataset_ids:
            for dataset_id in dataset_ids:
                if not validations.is_valid_dataset_id(dataset_id):
                    raise ValueError(f"Invalid dataset ID: {dataset_id}")
        if reaction_ids:
            for reaction_id in reaction_ids:
                if not validations.is_valid_reaction_id(reaction_id):
                    raise ValueError(f"Invalid reaction ID: {reaction_id}")
        if components is not None:
            component_params = [component.get_params() for component in components]
        else:
            component_params = None
        params = {
            "dataset_ids": ",".join(dataset_ids) if dataset_ids else None,
            "reaction_ids": ",".join(reaction_ids) if reaction_ids else None,
            "reaction_smarts": reaction_smarts,
            "dois": ",".join(dois) if dois else None,
            "component": component_params,
            "use_stereochemistry": use_stereochemistry,
            "similarity": similarity,
        }
        target = self._target + self._prefix + "/api/query"
        response = requests.get(target, params=params, timeout=timeout)
        return parse_response(response)


class ComponentQuery:
    """Client-side implementation of ReactionComponentPredicate."""

    def __init__(self, pattern: str, target: str, mode: str):
        """Initializes the query.

        Args:
            pattern: String SMILES or SMARTS pattern.
            target: String ReactionComponentPredicate.Target value.
            mode: String ReactionComponentPredicate.MatchMode value.
        """
        self._pattern = pattern
        try:
            query.ReactionComponentPredicate.Target.from_name(target)
        except KeyError as error:
            raise ValueError(f"target is not in {query.ReactionComponentPredicate.Target}: {target}") from error
        self._target = target
        try:
            query.ReactionComponentPredicate.MatchMode.from_name(mode)
        except KeyError as error:
            raise ValueError(f"mode is not in {query.ReactionComponentPredicate.MatchMode}: {mode}") from error
        self._mode = mode

    def get_params(self):
        """Returns URL parameters for GET requests."""
        return ";".join([self._pattern, self._target, self._mode])
