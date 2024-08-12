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

"""Celery configuration."""

import os

HOST = os.environ.get("REDIS_HOST", "localhost")
PORT = os.environ.get("REDIS_PORT", "6379")

# pylint: disable=invalid-name

broker_url = f"redis://{HOST}:{PORT}/1"
imports = ("ord_interface.api.search",)
result_backend = f"redis://{HOST}:{PORT}/2"
