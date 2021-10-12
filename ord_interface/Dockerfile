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
# $ docker build \
#     --file=Dockerfile \
#     -t openreactiondatabase/ord-interface \
#     ..
#
# To push the built image to Docker Hub:
# docker push openreactiondatabase/ord-interface

FROM python:3.8

RUN apt-get update \
 && apt-get install -y \
    build-essential \
    libpq-dev \
    procps \
    python3-dev \
    unzip \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade pip \
 && pip install gunicorn

# Install ord-schema.
WORKDIR /usr/src/app
RUN git clone https://github.com/Open-Reaction-Database/ord-schema.git
WORKDIR ord-schema
ARG ORD_SCHEMA_TAG=v0.3.13
RUN git fetch --tags && git checkout "${ORD_SCHEMA_TAG}"
RUN pip install -r requirements.txt
RUN python setup.py install

# COPY the local state.
WORKDIR ../ord-interface
COPY requirements.txt setup.py ./
COPY ord_interface/ ord_interface/
RUN pip install -r requirements.txt
RUN python setup.py install

# Fix some missing deps in ord-schema.
# TODO(skearnes): Move into ord-schema requirements.txt.
RUN pip install pillow rdkit-pypi

# Build and launch the interface.
EXPOSE 5000
CMD gunicorn ord_interface.search:app \
    --bind 0.0.0.0:5000 \
    --workers 2 \
    --access-logfile - \
    --access-logformat '%(t)s %({user-id}o)s %(U)s %(s)s %(L)s %(b)s %(f)s "%(r)s" "%(a)s"'
