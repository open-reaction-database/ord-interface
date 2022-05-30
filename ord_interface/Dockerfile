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

# To build the interface using the current state of ord-interface:
# $ cd path/to/ord-interface/ord_interface
# $ docker build --file=Dockerfile -t openreactiondatabase/ord-interface ..

FROM continuumio/miniconda3

# default-jre is required for running the closure compiler linter.
# For this next line, see:
# https://github.com/geerlingguy/ansible-role-java/issues/64#issuecomment-597132394
RUN mkdir -p /usr/share/man/man1/
RUN apt-get update \
 && apt-get install -y \
    build-essential \
    default-jre \
    npm \
    procps \
    unzip \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN conda install -c rdkit \
    psycopg2 \
    python=3.7 \
    rdkit \
 && conda clean -afy

# Fetch and build editor dependencies.
# NOTE(kearnes): Do this before COPYing the local state so it can be cached.
WORKDIR /usr/src/app/ord-interface/ord_interface/editor
# Bust the cache to get the latest commits; see https://stackoverflow.com/a/39278224.
ADD "https://api.github.com/repos/Open-Reaction-Database/ketcher/git/refs/heads/main" ketcher-version.json
RUN git clone https://github.com/Open-Reaction-Database/ketcher.git
RUN cd ketcher \
 && rm package-lock.json \
 && npm install \
 && npm run build \
 && rm -rf node_modules
RUN wget https://github.com/google/closure-library/archive/v20200517.tar.gz \
 && tar -xzf v20200517.tar.gz \
 && rm v20200517.tar.gz
RUN wget https://github.com/protocolbuffers/protobuf/releases/download/v3.14.0/protobuf-js-3.14.0.tar.gz \
 && tar -xzf protobuf-js-3.14.0.tar.gz \
 && rm protobuf-js-3.14.0.tar.gz
RUN wget https://raw.githubusercontent.com/google/closure-compiler/master/contrib/externs/jquery-3.3.js \
 && mkdir -p externs \
 && mv jquery-3.3.js externs
RUN npm install google-closure-compiler

# Install ord-schema.
WORKDIR /usr/src/app
RUN git clone https://github.com/Open-Reaction-Database/ord-schema.git
WORKDIR ord-schema
ARG ORD_SCHEMA_TAG=v0.3.15
RUN git fetch --tags && git checkout "${ORD_SCHEMA_TAG}"
RUN pip install -r requirements.txt
RUN python setup.py install

# Install interface dependencies.
WORKDIR ../ord-interface
COPY requirements.txt ./
RUN pip install -r requirements.txt
RUN pip install gunicorn

# Copy the local state.
COPY ord_interface/editor/js/ ord_interface/editor/js/
COPY ord_interface/editor/Makefile ord_interface/editor
WORKDIR ord_interface/editor
RUN make
WORKDIR ../..
# Copy the rest of the editor files.
COPY ord_interface/editor/ ord_interface/editor/
# Now copy everything else.
COPY ord_interface/interface.py ord_interface/
COPY ord_interface/client/ ord_interface/client/
COPY ord_interface/visualization/ ord_interface/visualization/
COPY setup.py .
RUN python setup.py install

# Build and launch the interface.
EXPOSE 5000
CMD gunicorn ord_interface.interface:app \
    --bind 0.0.0.0:5000 \
    --workers 2 \
    --access-logfile - \
    --access-logformat '%(t)s %({user-id}o)s %(U)s %(s)s %(L)s %(b)s %(f)s "%(r)s" "%(a)s"'
