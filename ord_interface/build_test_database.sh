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

export PGPASSWORD=postgres
# Use a non-standard PGDATA so the database persists; see
# https://nickjanetakis.com/blog/docker-tip-79-saving-a-postgres-database-in-a-docker-image.
CONTAINER="$(docker run --rm -d -p 5432:5432 -e POSTGRES_PASSWORD=${PGPASSWORD} -e PGDATA=/data informaticsmatters/rdkit-cartridge-debian)"

# Wait for the database to become available.
function connect() {
  psql -p 5432 -h localhost -U postgres < /dev/null
}
set +e
connect
while [ $? -ne 0 ]; do
  echo waiting for postgres
  sleep 1
  connect
done
set -e

# Editor.
cd editor
psql -p 5432 -h localhost -U postgres -f schema.sql
python py/migrate.py
cd ..

# Client.
psql -p 5432 -h localhost -U postgres -c 'CREATE DATABASE ord;'
wget --no-clobber https://github.com/open-reaction-database/ord-data/raw/main/data/89/ord_dataset-89b083710e2d441aa0040c361d63359f.pb.gz
wget --no-clobber https://github.com/open-reaction-database/ord-data/raw/main/data/b4/ord_dataset-b440f8c90b6343189093770060fc4098.pb.gz
python client/build_database.py --input="*.pb.gz"

# Save and shut down the container.
docker commit "${CONTAINER}" "openreactiondatabase/ord-postgres:test"
docker stop "${CONTAINER}"
