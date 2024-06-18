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
"""Populates a PostgreSQL database containing the ORD (for local testing)."""
import logging

from ord_schema.orm import database

from ord_interface.api.testing import setup_test_postgres
from ord_interface.client import POSTGRES_PORT, POSTGRES_PASSWORD, POSTGRES_DB, POSTGRES_USER, POSTGRES_HOST


def main():
    logging.disable(logging.DEBUG)
    connection_string = database.get_connection_string(
        database=POSTGRES_DB, username=POSTGRES_USER, password=POSTGRES_PASSWORD, host=POSTGRES_HOST, port=POSTGRES_PORT
    )
    setup_test_postgres(connection_string)


if __name__ == "__main__":
    main()
