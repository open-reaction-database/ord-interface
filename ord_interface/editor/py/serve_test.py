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
"""Tests for editor.py.serve."""
import base64
import json
import os
from urllib import parse

import flask
import pytest
from google.protobuf import text_format
from rdkit import Chem

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

import serve  # pylint: disable=import-error,wrong-import-order

# These temporary datasets are leaked by tests and must be deleted in setUp().
DATASETS = [
    "dataset",
    "other",
    "test",
]
TESTDATA = os.path.join(os.path.dirname(os.path.realpath(__file__)), "testdata")


@pytest.fixture
def client():
    app = flask.Flask(__name__)
    app.register_blueprint(serve.bp, url_prefix="")
    app.testing = True
    client = app.test_client()
    # GET requests automatically login as the test user.
    client.get("/authenticate")
    _destroy_datasets(client)
    # Start with an initial empty dataset called 'dataset'.
    _create(client, "dataset")
    yield client


def _create(client, dataset) -> None:
    """Make an empty dataset for testing."""
    response = client.post(f"/dataset/{dataset}/new", follow_redirects=True)
    assert response.status_code == 200


def _destroy(client, dataset) -> None:
    """Clean up a dataset that was created for testing."""
    response = client.get(f"/dataset/{dataset}/delete", follow_redirects=True)
    # Returns 200 even if the dataset did not exist before.
    assert response.status_code == 200


def _destroy_datasets(client) -> None:
    for dataset in DATASETS:
        _destroy(client, dataset)


def _get_dataset() -> dataset_pb2.Dataset:
    """Returns a Dataset for testing."""
    dataset = dataset_pb2.Dataset()
    with open(os.path.join(TESTDATA, "nielsen_fig1_dataset.pbtxt"), "rt") as f:
        text_format.Parse(f.read(), dataset)
    # Add some unicode to check for encoding/decoding robustness.
    # From https://en.wikipedia.org/wiki/Atlantis.
    dataset.reactions[0].provenance.city = "Ἀτλαντὶς νῆσος"
    return dataset


def _download_dataset(client, name) -> dataset_pb2.Dataset:
    """Downloads an existing dataset."""
    response = client.get(f"/dataset/{name}/download", follow_redirects=True)
    assert response.status_code == 200
    return dataset_pb2.Dataset.FromString(response.data)


def _upload_dataset(client, dataset, name) -> None:
    """Uploads a Dataset for testing."""
    response = client.post(f"/dataset/{name}/upload", data=text_format.MessageToString(dataset), follow_redirects=True)
    assert response.status_code == 200


def test_show_root(client):
    response = client.get("/", follow_redirects=True)
    assert response.status_code == 200


def test_show_datasets(client):
    response = client.get("/datasets", follow_redirects=True)
    assert response.status_code == 200


@pytest.mark.parametrize("dataset,expected", (("dataset", 200), ("other", 404)))
def test_show_dataset(client, dataset, expected):
    response = client.get(f"/dataset/{dataset}", follow_redirects=True)
    assert response.status_code == expected


@pytest.mark.parametrize(
    "filename,expected",
    (
        ("dataset", 200),
        ("../dataset", 404),
        (parse.quote_plus("../dataset"), 404),
        ("/foo/bar", 404),
        ("other", 404),
    ),
)
def test_download_dataset(client, filename, expected, tmp_path):
    response = client.get(f"/dataset/{filename}/download", follow_redirects=True)
    assert response.status_code == expected
    if response.status_code == 200:
        # Make sure it parses.
        filename = (tmp_path / "dataset.pb").as_posix()
        with open(filename, "wb") as f:
            f.write(response.data)
        message_helpers.load_message(filename, dataset_pb2.Dataset)


@pytest.mark.parametrize(
    "filename,kind,expected",
    (
        ("dataset", "pb", 200),
        ("dataset", "pbtxt", 200),
        ("../dataset", "pb", 404),
        ("../dataset", "pbtxt", 404),
        (parse.quote_plus("../dataset"), "pb", 404),
        (parse.quote_plus("../dataset"), "pbtxt", 404),
        ("/foo/bar", "pb", 404),
        ("/foo/bar", "pbtxt", 404),
        ("other", "pb", 404),
        ("other", "pbtxt", 404),
    ),
)
def test_download_dataset_with_kind(client, filename, kind, expected, tmp_path):
    response = client.get(f"/dataset/{filename}/download/{kind}", follow_redirects=True)
    assert response.status_code == expected
    if response.status_code == 200:
        # Make sure it parses.
        filename = (tmp_path / f"dataset.{kind}").as_posix()
        with open(filename, "wb") as f:
            f.write(response.data)
        message_helpers.load_message(filename, dataset_pb2.Dataset)


@pytest.mark.parametrize(
    "filename,expected,as_text",
    (
        ("dataset", 409, True),
        ("dataset", 409, False),
        ("other", 200, True),
        ("other", 200, False),
    ),
)
def test_upload_dataset(client, filename, expected, as_text):
    dataset = _get_dataset()
    if as_text:
        data = text_format.MessageToString(dataset)
    else:
        data = dataset.SerializeToString()
    response = client.post(f"/dataset/{filename}/upload", data=data, follow_redirects=True)
    assert response.status_code == expected
    if response.status_code == 200:
        response = client.get(f"/dataset/{filename}/download", follow_redirects=True)
        assert response.status_code == 200
        downloaded_dataset = dataset_pb2.Dataset.FromString(response.data)
        assert downloaded_dataset == dataset


@pytest.mark.parametrize("name,expected", (("dataset", 409), ("other", 200)))
def test_new_dataset(client, name, expected):
    response = client.post(f"/dataset/{name}/new", follow_redirects=True)
    assert response.status_code == expected
    if response.status_code == 200:
        dataset = _download_dataset(client, name)
        assert len(dataset.reactions) == 0


@pytest.mark.parametrize("prefix", (b"", b"data:foo/bar;base64,"))
def test_enumerate_dataset(client, prefix):
    data = {"spreadsheet_name": "test.csv"}
    with open(os.path.join(TESTDATA, "nielsen_fig1.csv"), "rb") as f:
        data["spreadsheet_data"] = (prefix + base64.b64encode(f.read())).decode()
    with open(os.path.join(TESTDATA, "nielsen_fig1_template.pbtxt"), "rt") as f:
        data["template_string"] = f.read()
    response = client.post("/dataset/enumerate", json=data, follow_redirects=True)
    assert response.status_code == 200
    response = client.get("/dataset/test_dataset/download", follow_redirects=True)
    assert response.status_code == 200
    dataset = dataset_pb2.Dataset.FromString(response.data)
    assert len(dataset.reactions) == 80


@pytest.mark.parametrize("index,expected", ((0, 200), (3, 200), (80, 404)))
def test_show_reaction(client, index, expected):
    _upload_dataset(client, _get_dataset(), "test")
    response = client.get(f"/dataset/test/reaction/{index}", follow_redirects=True)
    assert response.status_code == expected


def test_download_reaction(client):
    reaction = _get_dataset().reactions[0]
    response = client.post("/reaction/download", data=reaction.SerializeToString(), follow_redirects=True)
    assert response.status_code == 200
    downloaded_reaction = reaction_pb2.Reaction()
    text_format.Parse(response.data, downloaded_reaction)
    assert downloaded_reaction == reaction


def test_new_reaction(client):
    name = "test"
    dataset = _get_dataset()
    _upload_dataset(client, dataset, name)
    response = client.get(f"/dataset/{name}/new/reaction", follow_redirects=True)
    assert response.status_code == 200
    downloaded_dataset = _download_dataset(client, name)
    assert len(downloaded_dataset.reactions) == 81


def test_clone_reaction(client):
    name = "test"
    dataset = _get_dataset()
    _upload_dataset(client, dataset, name)
    response = client.get(f"/dataset/{name}/clone/0", follow_redirects=True)
    assert response.status_code == 200
    downloaded_dataset = _download_dataset(client, name)
    assert len(downloaded_dataset.reactions) == 81
    assert dataset.reactions[0] == downloaded_dataset.reactions[80]


def test_delete_reaction(client):
    name = "test"
    dataset = _get_dataset()
    _upload_dataset(client, dataset, name)
    response = client.get(f"/dataset/{name}/delete/reaction/0", follow_redirects=True)
    assert response.status_code == 200
    downloaded_dataset = _download_dataset(client, name)
    assert len(downloaded_dataset.reactions) == 79
    assert dataset.reactions[1] == downloaded_dataset.reactions[0]


def test_delete_reaction_id(client):
    name = "test"
    dataset = dataset_pb2.Dataset()
    reaction_id = "test_reaction_id"
    dataset.reaction_ids.append(reaction_id)
    _upload_dataset(client, dataset, name)
    response = client.get(f"/dataset/{name}/delete/reaction_id/{reaction_id}", follow_redirects=True)
    assert response.status_code == 200
    downloaded_dataset = _download_dataset(client, name)
    assert len(downloaded_dataset.reaction_ids) == 0


def test_delete_reaction_id_blank(client):
    name = "test"
    dataset = dataset_pb2.Dataset(reaction_ids=["", "test", ""])
    _upload_dataset(client, dataset, name)
    response = client.get(f"/dataset/{name}/delete/reaction_id", follow_redirects=True)
    assert response.status_code == 200
    downloaded_dataset = _download_dataset(client, name)
    assert len(downloaded_dataset.reaction_ids) == 2


def test_read_dataset(client):
    name = "test"
    dataset = _get_dataset()
    _upload_dataset(client, dataset, name)
    response = client.get(f"/dataset/proto/read/{name}", follow_redirects=True)
    assert response.status_code == 200
    downloaded_dataset = dataset_pb2.Dataset()
    downloaded_dataset.ParseFromString(response.data)
    assert downloaded_dataset == dataset


def test_write_dataset(client):
    name = "test"
    dataset = _get_dataset()
    response = client.post(f"/dataset/proto/write/{name}", data=dataset.SerializeToString(), follow_redirects=True)
    assert response.status_code == 200
    downloaded_dataset = _download_dataset(client, name)
    assert downloaded_dataset == dataset


def test_write_upload(client):
    name = "test"
    data = b"test data"
    token = b"upload_token"
    dataset = dataset_pb2.Dataset()
    reaction = dataset.reactions.add()
    observation = reaction.observations.add()
    observation.image.bytes_value = token
    _upload_dataset(client, dataset, name)
    response = client.post(f"/dataset/proto/upload/{name}/{token.decode()}", data=data, follow_redirects=True)
    assert response.status_code == 200
    # Verify that the token was resolved in the Dataset.
    downloaded_dataset = _download_dataset(client, name)
    assert downloaded_dataset.reactions[0].observations[0].image.bytes_value == data


def test_read_upload(client):
    data = b"test data"
    token = "upload_token"
    response = client.post(f"/dataset/proto/download/{token}", data=data, follow_redirects=True)
    assert response.status_code == 200
    assert response.data == data


@pytest.mark.parametrize(
    "message,expected_num_errors,expected_num_warnings",
    (
        (reaction_pb2.Percentage(value=15.6), 0, 0),
        (reaction_pb2.Percentage(precision=-15.6), 2, 0),
    ),
)
def test_validate_reaction(client, message, expected_num_errors, expected_num_warnings):
    response = client.post(
        f"/dataset/proto/validate/{message.DESCRIPTOR.name}", data=message.SerializeToString(), follow_redirects=True
    )
    assert response.status_code == 200
    output = json.loads(response.data)
    assert len(output["errors"]) == expected_num_errors
    assert len(output["warnings"]) == expected_num_warnings


@pytest.mark.parametrize("identifier_type,data,expected", (("NAME", "benzene", "c1ccccc1"),))
def test_resolve_compound(client, identifier_type, data, expected):
    response = client.post(f"/resolve/{identifier_type}", data=data, follow_redirects=True)
    assert response.status_code == 200
    resolved, _ = json.loads(response.data)
    # NOTE(kearnes): Try to compensate for values from different services.
    canonical_resolved = Chem.MolToSmiles(Chem.MolFromSmiles(resolved))
    assert canonical_resolved == expected


def test_render_reaction(client):
    reaction = reaction_pb2.Reaction()
    component = reaction.inputs["test"].components.add()
    component.identifiers.add(value="c1ccccc1", type="SMILES")
    response = client.post("/render/reaction", data=reaction.SerializeToString(), follow_redirects=True)
    assert response.status_code == 200


def test_render_compound(client):
    compound = reaction_pb2.Compound()
    compound.identifiers.add(value="c1ccccc1", type="SMILES")
    response = client.post("/render/reaction", data=compound.SerializeToString(), follow_redirects=True)
    assert response.status_code == 200


def test_compare(client):
    name = "test"
    dataset = _get_dataset()
    _upload_dataset(client, dataset, name)
    response = client.post(f"/dataset/proto/compare/{name}", data=dataset.SerializeToString(), follow_redirects=True)
    assert response.status_code == 200
    dataset.reactions[0].reaction_id = "not the original"
    response = client.post(f"/dataset/proto/compare/{name}", data=dataset.SerializeToString(), follow_redirects=True)
    assert response.status_code == 409


@pytest.mark.skip("Requires the editor to be built.")
def test_js():
    pass


@pytest.mark.parametrize("sheet,expected", (("reaction.css", 200), ("percentage.css", 404)))
def test_css(client, sheet, expected):
    response = client.get(f"/css/{sheet}", follow_redirects=True)
    assert response.status_code == expected


@pytest.mark.parametrize("path", ("dataset/deps.js", "dataset/test/deps.js", "dataset/test/reaction/deps.js"))
def test_deps(client, path):
    response = client.get(path, follow_redirects=True)
    assert response.status_code == 200


def test_get_molfile(client):
    smiles = "c1ccccc1"
    compound = reaction_pb2.Compound()
    compound.identifiers.add(value=smiles, type="SMILES")
    response = client.post("/ketcher/molfile", data=compound.SerializeToString(), follow_redirects=True)
    assert response.status_code == 200
    assert json.loads(response.data) == Chem.MolToMolBlock(Chem.MolFromSmiles(smiles))


def test_get_molfile_no_structure(client):
    compound = reaction_pb2.Compound()
    compound.identifiers.add(value="benzene", type="NAME")
    response = client.post("/ketcher/molfile", data=compound.SerializeToString(), follow_redirects=True)
    assert response.status_code == 204
