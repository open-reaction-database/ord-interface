#!/bin/bash
# Copyright 2022 Open Reaction Database Project Authors
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

set -e

set -x
./build_test_database.sh
ARCH="x86_64"
if [[ "$(uname -p)" =~ "arm" ]]; then
  ARCH="aarch_64"
fi
docker build -f Dockerfile -t openreactiondatabase/ord-interface .. --build-arg="ARCH=${ARCH}" "$@"
docker compose up --detach
set +x

# Wait for the database to become available.
function connect() {
  PGPASSWORD=postgres psql -p 5432 -h localhost -U postgres < /dev/null
}
set +e
connect
while [ $? -ne 0 ]; do
  echo waiting for postgres
  sleep 1
  connect
done
status=0

# Run tests (only the ones that depend on the running app; not those that use testing.postgresql).
export POSTGRES_HOST=localhost
export POSTGRES_USER=postgres
export POSTGRES_PASSWORD=postgres
export POSTGRES_DATABASE=ord
pytest -vv --ignore=api/queries_test.py || status=1
node editor/js/test.js || status=1

# Shut down the containers.
docker compose down

# Report pass/fail.
red='\033[0;31m'
green='\033[0;32m'
neutral='\033[0m'
[ "${status}" -eq 0 ] && \
    printf "${green}PASS${neutral}\n" || printf "${red}FAIL${neutral}\n"

# Relay the status for GitHub CI.
test "${status}" -eq 0
