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
"""Query web interface to the Open Reaction Database in Postgres.

The client is stateless. Full information is in the URL so results can be
linked. Query parameters are communicated in URL GET params:

  component=<pattern;source;(exact|substructure|similarity|smarts)>

    The second token specifies whether the predicate should match an input or
    an output.

    The last token specifies the matching criterion. The default is "exact".
    The pattern is a SMILES string, unless the token is "smarts" in which case
    the pattern is a SMARTS string.

    Component may be repeated any number of times.

  reaction_ids=<ids>

  reaction_smarts=<smarts>

These query types are mutually exclusive. All query parameters are assumed to
be URL-encoded.
"""

# pylint: disable=too-many-locals

import dataclasses
import os
from typing import Dict, List, Optional, Tuple, Union

import flask
from rdkit import Chem

from ord_interface.client import query
from ord_interface.visualization import generate_text

bp = flask.Blueprint("client", __name__, url_prefix="/client", template_folder=".")
POSTGRES_HOST = os.getenv("POSTGRES_HOST", "localhost")
POSTGRES_PORT = os.getenv("POSTGRES_PORT", "5432")
POSTGRES_USER = os.getenv("POSTGRES_USER", "ord-postgres")
POSTGRES_PASSWORD = os.getenv("POSTGRES_PASSWORD", "ord-postgres")
POSTGRES_DATABASE = os.getenv("POSTGRES_DATABASE", "ord")

BOND_LENGTH = 20


@bp.route("/")
def show_root():
    flask.redirect(flask.url_for(".show_browse"))


@bp.route("/browse")
def show_browse():
    """Shows the browser interface."""
    return flask.render_template("browse.html", datasets=fetch_datasets())


@bp.route("/search")
def show_search():
    """Shows the search interface.

    Creates a query to show a set of randomly selected reactions so the
    page won't be empty.
    """
    command, limit = build_query()
    if command is None:
        command = query.RandomSampleQuery(100)
    query_json = command.json()
    try:
        results = connect().run_query(command, limit=limit, return_ids=True)
        error = None
    except query.QueryException as exception:
        results = None
        error = f"(Error) {exception}"
    if results is not None and not results:
        results = None
        error = "query did not match any reactions"
    return flask.render_template("search.html", results=results, error=error, query=query_json)


@bp.route("/id/<reaction_id>")
def show_id(reaction_id):
    """Returns the pbtxt of a single reaction as plain text."""
    results = connect().run_query(query.ReactionIdQuery([reaction_id]))
    if len(results) == 0:
        return flask.abort(404)
    reaction = results[0].reaction
    reaction_summary = generate_text.generate_html(reaction, bond_length=BOND_LENGTH)
    return flask.render_template(
        "reaction_view.html",
        reaction=reaction,
        dataset_id=results[0].dataset_id,
        reaction_summary=reaction_summary,
        bond_length=BOND_LENGTH,
    )


@bp.route("/render/<reaction_id>")
def render_reaction(reaction_id):
    """Renders a reaction as an HTML table with images and text."""
    command = query.ReactionIdQuery([reaction_id])
    results = connect().run_query(command)
    if len(results) == 0 or len(results) > 1:
        return flask.abort(404)
    result = results[0]
    try:
        html = generate_text.generate_html(reaction=result.reaction, compact=True)
        return flask.jsonify(html)
    except (ValueError, KeyError):
        return flask.jsonify("[Reaction cannot be displayed]")


def connect():
    return query.OrdPostgres(
        dbname=POSTGRES_DATABASE,
        user=POSTGRES_USER,
        password=POSTGRES_PASSWORD,
        host=POSTGRES_HOST,
        port=int(POSTGRES_PORT),
    )


@bp.route("/api/fetch_reactions", methods=["POST"])
def fetch_reactions():
    """Fetches a list of Reactions by ID."""
    reaction_ids = flask.request.get_json()
    command = query.ReactionIdQuery(reaction_ids)
    try:
        results = connect().run_query(command)
        return flask.jsonify([dataclasses.asdict(result) for result in results])
    except query.QueryException as error:
        return flask.abort(flask.make_response(str(error), 400))


def fetch_datasets() -> List[Dict[str, Union[str, int]]]:
    """Fetches info about the current datasets."""
    engine = connect()
    rows = {}
    with engine.connection, engine.cursor() as cursor:
        cursor.execute("SELECT dataset_id, name, description FROM datasets")
        for dataset_id, name, description in cursor:
            rows[dataset_id] = {
                "Dataset ID": dataset_id,
                "Name": name,
                "Description": description,
                "Size": 0,
            }
        # Get dataset sizes.
        cursor.execute(
            """
            SELECT dataset_id, COUNT(reaction_id)
            FROM reactions
            GROUP BY dataset_id
            """
        )
        for dataset_id, count in cursor:
            rows[dataset_id]["Size"] = count
        return list(rows.values())


@bp.route("/api/query")
def run_query():
    """Builds and executes a GET query.

    Returns:
        A serialized Dataset proto containing the matched reactions.
    """
    command, limit = build_query()
    if command is None:
        return flask.abort(flask.make_response("no query defined", 400))
    try:
        results = connect().run_query(command, limit=limit)
        return flask.jsonify([dataclasses.asdict(result) for result in results])
    except query.QueryException as error:
        return flask.abort(flask.make_response(str(error), 400))


def build_query() -> Tuple[Optional[query.ReactionQueryBase], Optional[int]]:
    """Builds a query from GET parameters.

    Returns:
        query: ReactionQueryBase subclass instance.
        limit: Maximum number of results to return.
    """
    dataset_ids = flask.request.args.get("dataset_ids")
    reaction_ids = flask.request.args.get("reaction_ids")
    reaction_smarts = flask.request.args.get("reaction_smarts")
    dois = flask.request.args.get("dois")
    components = flask.request.args.getlist("component")
    use_stereochemistry = flask.request.args.get("use_stereochemistry")
    similarity = flask.request.args.get("similarity")
    limit = flask.request.args.get("limit")
    if limit is not None:
        limit = int(limit)
    if dataset_ids is not None:
        command = query.DatasetIdQuery(dataset_ids.split(","))
    elif reaction_ids is not None:
        command = query.ReactionIdQuery(reaction_ids.split(","))
    elif reaction_smarts is not None:
        command = query.ReactionSmartsQuery(reaction_smarts)
    elif dois is not None:
        command = query.DoiQuery(dois.split(","))
    elif components:
        predicates = []
        for component in components:
            pattern, source, mode_name = component.split(";")
            table = query.ReactionComponentPredicate.SOURCE_TO_TABLE[source]
            mode = query.ReactionComponentPredicate.MatchMode.from_name(mode_name)
            predicates.append(query.ReactionComponentPredicate(pattern, table, mode))
        kwargs = {}
        if use_stereochemistry is not None:
            kwargs["do_chiral_sss"] = use_stereochemistry
        if similarity is not None:
            kwargs["tanimoto_threshold"] = float(similarity)
        command = query.ReactionComponentQuery(predicates, **kwargs)
    else:
        command = None
    return command, limit


@bp.route("/ketcher/molfile", methods=["POST"])
def get_molfile():
    """Returns a molblock for the given SMILES."""
    smiles = flask.request.get_data()
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(smiles)
        return flask.jsonify(Chem.MolToMolBlock(mol))
    except ValueError:
        return f"could not parse SMILES: {smiles}", 400
