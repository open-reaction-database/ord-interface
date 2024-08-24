# Copyright 2024 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Editor API helper functions."""
import gzip
from base64 import b64encode

from fastapi import Response
from google.protobuf import json_format, text_format  # pytype: disable=import-error
from ord_schema.proto.dataset_pb2 import Dataset
from ord_schema.proto.reaction_pb2 import Reaction

from ord_interface.editor.api.database import add_dataset, get_cursor, get_dataset
from ord_interface.editor.api.testing import setup_test_postgres


def download_message(message: Dataset | Reaction, prefix: str, kind: str):
    """Downloads a dataset or reaction as a gzipped file."""
    match kind:
        case "json":
            data = json_format.MessageToJson(message)
            filename = f"{prefix}.json"
        case "binpb":
            data = message.SerializeToString()
            filename = f"{prefix}.binpb"
        case "txtpb":
            data = text_format.MessageToBytes(message)
            filename = f"{prefix}.txtpb"
        case _:
            raise ValueError(kind)
    headers = {"Content-Disposition": f'attachment; filename="{filename}.gz"'}
    return Response(gzip.compress(data), headers=headers, media_type="application/gzip")


def send_message(message) -> str:
    """Converts a protocol buffer message to a base64-encoded string."""
    return b64encode(message.SerializeToString()).decode()
