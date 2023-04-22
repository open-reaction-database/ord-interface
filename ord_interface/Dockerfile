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

FROM python:3.10

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

# Fetch and build editor dependencies.
# NOTE(kearnes): Do this before COPYing the local state so it can be cached.
WORKDIR /usr/src/app/ord-interface/ord_interface
RUN wget https://github.com/epam/ketcher/releases/download/v2.5.1/ketcher-standalone-2.5.1.zip \
 && unzip ketcher-standalone-2.5.1.zip
WORKDIR editor
RUN wget https://github.com/google/closure-library/archive/v20200517.tar.gz \
 && tar -xzf v20200517.tar.gz \
 && rm v20200517.tar.gz
RUN wget https://github.com/protocolbuffers/protobuf/releases/download/v3.14.0/protobuf-js-3.14.0.tar.gz \
 && tar -xzf protobuf-js-3.14.0.tar.gz \
 && rm protobuf-js-3.14.0.tar.gz
RUN wget https://github.com/protocolbuffers/protobuf/releases/download/v3.14.0/protoc-3.14.0-linux-x86_64.zip \
 && unzip protoc-3.14.0-linux-x86_64.zip
ENV PATH="/usr/src/app/ord-interface/ord_interface/editor/bin:${PATH}"
RUN wget https://raw.githubusercontent.com/google/closure-compiler/master/contrib/externs/jquery-3.3.js \
 && mkdir -p externs \
 && mv jquery-3.3.js externs
ARG ORD_SCHEMA=0.3.49
RUN wget https://github.com/open-reaction-database/ord-schema/archive/refs/tags/v${ORD_SCHEMA}.tar.gz \
 && tar -xzf v${ORD_SCHEMA}.tar.gz \
 && mv ord-schema-${ORD_SCHEMA} ord-schema \
 && rm v${ORD_SCHEMA}.tar.gz
RUN npm install google-closure-compiler@20220803

# Install ord-interface.
WORKDIR /usr/src/app/ord-interface
COPY . .
RUN pip install --upgrade pip \
 && pip install gunicorn wheel \
 && pip install . \
 && pip cache purge
RUN cd ord_interface/editor \
 && ./build.sh \
 && cd ../..

# Build and launch the interface.
EXPOSE 5001
CMD gunicorn ord_interface.interface:app \
    --bind 0.0.0.0:5001 \
    --workers 2 \
    --access-logfile - \
    --access-logformat '%(t)s %({user-id}o)s %(U)s %(s)s %(L)s %(b)s %(f)s "%(r)s" "%(a)s"'
