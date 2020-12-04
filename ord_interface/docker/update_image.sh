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

set -ex

docker build \
  --file=ord_interface/docker/Dockerfile \
  -t ord-postgres:empty \
  .
CONTAINER="$(docker run --rm -d ord-postgres:empty)"
echo "Waiting 5s for the server to start..."
sleep 5
TAG="${ORD_POSTGRES_TAG:-latest}"
# NOTE(kearnes): Pass the tag directly instead of using ORD_POSTGRES_TAG
# so we don't have to set that variable in the docker image.
docker exec "${CONTAINER}" ./ord_interface/docker/build_database.sh "${TAG}"
docker commit "${CONTAINER}" "openreactiondatabase/ord-postgres:${TAG}"
docker stop "${CONTAINER}"
