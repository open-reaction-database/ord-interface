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

"""Open Reaction Database API."""

import os
from contextlib import asynccontextmanager

from fastapi import FastAPI
from testing.postgresql import Postgresql

from ord_interface.api import editor, search
from ord_interface.api.testing import setup_test_postgres


@asynccontextmanager
async def lifespan(*args, **kwargs):
    """FastAPI lifespan setup; see https://fastapi.tiangolo.com/advanced/events/#lifespan."""
    del args, kwargs  # Unused.
    if os.getenv("ORD_INTERFACE_TESTING", "FALSE") == "TRUE":
        with Postgresql() as postgres:
            setup_test_postgres(postgres.url())
            os.environ["ORD_INTERFACE_POSTGRES"] = postgres.url()
            yield
    else:
        yield


app = FastAPI(lifespan=lifespan)
app.include_router(search.router)
app.include_router(editor.router)


@app.get("/ketcher")
def ketcher():
    return app.redirect("/standalone/index.html")
