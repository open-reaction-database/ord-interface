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

import os
from typing import NewType, Optional, Tuple

import flask

from ord_schema.visualization import generate_text

from ord_interface import query

app = flask.Flask(__name__, template_folder='.')
app.config['ORD_POSTGRES_HOST'] = os.getenv('ORD_POSTGRES_HOST', 'localhost')

Query = NewType('Query', query.ReactionQueryBase)


@app.route('/')
def show_root():
    """Shows the web form.

    Creates a query to show a set of randomly selected reactions so the
    page won't be empty.
    """
    command, limit = build_query()
    if command is None:
        command = query.RandomSampleQuery(100)
    query_json = command.json()
    try:
        dataset = connect().run_query(command, limit=limit, return_ids=True)
        error = None
    except query.QueryException as exception:
        dataset = None
        error = f'(Error) {exception}'
    if dataset is not None and not dataset.reaction_ids:
        dataset = None
        error = 'query did not match any reactions'
    return flask.render_template('search.html',
                                 dataset=dataset,
                                 error=error,
                                 query=query_json)


@app.route('/id/<reaction_id>')
def show_id(reaction_id):
    """Returns the pbtxt of a single reaction as plain text."""
    dataset = connect().run_query(query.ReactionIdQuery([reaction_id]))
    if len(dataset.reactions) == 0:
        return flask.abort(404)
    return generate_text.generate_summary(dataset.reactions[0])


@app.route('/render/<reaction_id>')
def render_reaction(reaction_id):
    """Renders a reaction as an HTML table with images and text."""
    command = query.ReactionIdQuery([reaction_id])
    dataset = connect().run_query(command)
    if len(dataset.reactions) == 0 or len(dataset.reactions) > 1:
        return flask.abort(404)
    reaction = dataset.reactions[0]
    try:
        html = generate_text.generate_html(reaction, compact=True)
        return flask.jsonify(html)
    except (ValueError, KeyError):
        return flask.jsonify('[Reaction cannot be displayed]')


def connect():
    return query.OrdPostgres(dbname='ord',
                             user='ord-postgres',
                             password='ord-postgres',
                             host=app.config['ORD_POSTGRES_HOST'],
                             port=5432)


@app.route('/api/fetch_reactions', methods=['POST'])
def fetch_reactions():
    """Fetches a list of Reactions by ID."""
    reaction_ids = flask.request.get_json()
    command = query.ReactionIdQuery(reaction_ids)
    try:
        dataset = connect().run_query(command)
        return flask.make_response(dataset.SerializeToString())
    except query.QueryException as error:
        return flask.abort(flask.make_response(str(error), 400))


@app.route('/api/query')
def run_query():
    """Builds and executes a GET query.

    Returns:
        A serialized Dataset proto containing the matched reactions.
    """
    command, limit = build_query()
    if command is None:
        return flask.abort(flask.make_response('no query defined', 400))
    try:
        dataset = connect().run_query(command, limit=limit)
        return flask.make_response(dataset.SerializeToString())
    except query.QueryException as error:
        return flask.abort(flask.make_response(str(error), 400))


def build_query() -> Tuple[Optional[Query], Optional[int]]:
    """Builds a query from GET parameters.

    Returns:
        query: ReactionQueryBase subclass instance.
        limit: Maximum number of results to return.
    """
    dataset_ids = flask.request.args.get('dataset_ids')
    reaction_ids = flask.request.args.get('reaction_ids')
    reaction_smarts = flask.request.args.get('reaction_smarts')
    dois = flask.request.args.get('dois')
    components = flask.request.args.getlist('component')
    use_stereochemistry = flask.request.args.get('use_stereochemistry')
    similarity = flask.request.args.get('similarity')
    limit = flask.request.args.get('limit')
    if limit is not None:
        limit = int(limit)
    if dataset_ids is not None:
        command = query.DatasetIdQuery(dataset_ids.split(','))
    if reaction_ids is not None:
        command = query.ReactionIdQuery(reaction_ids.split(','))
    elif reaction_smarts is not None:
        command = query.ReactionSmartsQuery(reaction_smarts)
    elif dois is not None:
        command = query.DoiQuery(dois.split(','))
    elif components:
        predicates = []
        for component in components:
            pattern, source, mode_name = component.split(';')
            table = query.ReactionComponentPredicate.SOURCE_TO_TABLE[source]
            mode = query.ReactionComponentPredicate.MatchMode.from_name(
                mode_name)
            predicates.append(
                query.ReactionComponentPredicate(pattern, table, mode))
        kwargs = {}
        if use_stereochemistry is not None:
            kwargs['do_chiral_sss'] = use_stereochemistry
        if similarity is not None:
            kwargs['tanimoto_threshold'] = float(similarity)
        command = query.ReactionComponentQuery(predicates, **kwargs)
    else:
        command = None
    return command, limit
