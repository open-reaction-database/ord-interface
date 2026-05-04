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

FROM python:3.11

RUN mkdir -p /usr/share/man/man1/
RUN apt-get update \
   && apt-get install -y \
   build-essential \
   nginx \
   npm \
   procps \
   unzip \
   && apt-get clean \
   && rm -rf /var/lib/apt/lists/*

# Fetch the Ketcher distribution used by the Vue app.
# NOTE(kearnes): Do this before COPYing the local state so it can be cached.
WORKDIR /app/ord-interface/ord_interface
RUN wget https://github.com/epam/ketcher/releases/download/v2.5.1/ketcher-standalone-2.5.1.zip \
   && unzip ketcher-standalone-2.5.1.zip

# Install dependencies.
WORKDIR /app/ord-interface
COPY setup.py README.md ./
RUN pip install --upgrade pip \
   && pip install gunicorn wheel
# https://stackoverflow.com/a/60511098.
RUN python setup.py egg_info \
   && pip install $(grep -v '^\[' *.egg-info/requires.txt) \
   && rm -rf *.egg-info/

# Build the Vue app.
COPY app/ app/
WORKDIR /app/ord-interface/app
RUN cp -r /app/ord-interface/ord_interface/standalone/ /app/ord-interface/app/src/ketcher \
   && npm install \
   && npm run build \
   && mv dist/ /app/ord-interface/vue
WORKDIR /app/ord-interface
COPY ord_interface/api/ ord_interface/api/
COPY ord_interface/visualization/ ord_interface/visualization/
RUN pip install .

# Build and launch the interface.
COPY ord_interface/start_app.sh .
COPY ord_interface/nginx.conf /etc/nginx/nginx.conf

RUN ln -sf /dev/stderr /var/log/nginx.log \
   && chmod +x start_app.sh
ENTRYPOINT ["/app/ord-interface/start_app.sh"]
