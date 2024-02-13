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

# https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -Eeuxo pipefail

rm -rf gen  # Clean up from previous run.
mkdir -p gen/js/ord gen/js/proto/ord
# Compile JS proto wrappers.
# https://github.com/google/closure-library/releases/
CLOSURE=closure-library-20200517
# https://raw.githubusercontent.com/google/closure-compiler/master/contrib/externs/jquery-3.3.js
JQUERY_EXTERNS=externs/jquery-3.3.js
# https://github.com/protocolbuffers/protobuf/releases
PROTOBUF=protobuf-3.14.0
echo "protoc: $(which protoc) $(protoc --version)"
for source in /app/ord-schema/proto/*.proto; do
  protoc \
    --experimental_allow_proto3_optional \
    --js_out=binary:gen/js/proto/ord \
    --proto_path=/app \
    "${source}"
done
# Build dataset.js
java -jar node_modules/google-closure-compiler-java/compiler.jar \
  js/*.js \
  gen/js/proto/ord/*.js \
  "${PROTOBUF}"/js/*.js \
  "${PROTOBUF}"/js/binary/*.js \
  "${CLOSURE}"/closure/goog/*.js \
  "${CLOSURE}"/closure/goog/*/*.js \
  "${CLOSURE}"/closure/goog/labs/*/*.js \
  --externs="${JQUERY_EXTERNS}" \
  --entry_point=ord.dataset \
  --js_output_file=gen/js/ord/dataset.js \
  --dependency_mode=PRUNE \
  --jscomp_error="*" \
  --hide_warnings_for="${CLOSURE}" \
  --hide_warnings_for="${PROTOBUF}" \
  --hide_warnings_for=gen/js/proto/ord \
  || (rm -f gen/js/ord/dataset.js && false)
# Build reaction.js
java -jar node_modules/google-closure-compiler-java/compiler.jar \
  js/*.js \
  gen/js/proto/ord/*.js \
  "${PROTOBUF}"/js/*.js \
  "${PROTOBUF}"/js/binary/*.js \
  "${CLOSURE}"/closure/goog/*.js \
  "${CLOSURE}"/closure/goog/*/*.js \
  "${CLOSURE}"/closure/goog/labs/*/*.js \
  --externs="${JQUERY_EXTERNS}" \
  --entry_point=ord.reaction \
  --js_output_file=gen/js/ord/reaction.js \
  --dependency_mode=PRUNE \
  --jscomp_error="*" \
  --hide_warnings_for="${CLOSURE}" \
  --hide_warnings_for="${PROTOBUF}" \
  --hide_warnings_for=gen/js/proto/ord \
  || (rm -f gen/js/ord/reaction.js && false)
