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
"""Entrypoint for the web interface."""
import flask

from ord_interface.client import search
from ord_interface.editor.py import serve
from ord_interface.visualization import filters

app = flask.Flask(__name__)
# TODO(skearnes): Figure out bp.add_app_template_filter?
app.jinja_env.filters.update(filters.TEMPLATE_FILTERS)  # pylint: disable=no-member
app.register_blueprint(search.bp)
app.register_blueprint(serve.bp)


@app.route('/')
def show_root():
    return flask.redirect(search.bp.url_prefix)
