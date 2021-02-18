#!/bin/bash
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

set -e

# Build the ord-postgres image.
export ORD_POSTGRES_TAG=test
./ord_interface/docker/update_image.sh "$@"

# Build the ord-interface image and start the web server.
set -x
cd ord_interface
docker build -f Dockerfile -t openreactiondatabase/ord-interface ../ "$@"
docker-compose up --detach
set +x

echo "Waiting 5s for the server to start..."
sleep 5

set +e
status=0

CONTAINER="$(docker ps -q --filter name=web)"
docker exec "${CONTAINER}" python ord_interface/build_database_test.py || status=1
docker exec "${CONTAINER}" python ord_interface/ord_client_test.py || status=1
docker exec "${CONTAINER}" python ord_interface/query_test.py || status=1

# Shut down the containers.
docker-compose down

# Report pass/fail.
red='\033[0;31m'
green='\033[0;32m'
neutral='\033[0m'
[ "${status}" -eq 0 ] && \
    printf "${green}PASS${neutral}\n" || printf "${red}FAIL${neutral}\n"

# Relay the status for GitHub CI.
test "${status}" -eq 0
