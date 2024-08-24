# Copyright 2024 Open Reaction Database Project Authors
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

"""Tests for ord_interface.editor.api.database."""
from ord_interface.editor.api.database import get_dataset
from ord_interface.editor.api.testing import TEST_USER_ID


def test_get_dataset(test_cursor):
    dataset = get_dataset(TEST_USER_ID, "Deoxyfluorination screen", test_cursor)
    assert len(dataset.reactions) == 80


def test_get_unknown_dataset(test_cursor):
    dataset = get_dataset(TEST_USER_ID, "UNKNOWN", test_cursor)
    assert dataset is None
