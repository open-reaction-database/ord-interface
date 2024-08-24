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

"""User API endpoints."""

from uuid import uuid4

from fastapi import APIRouter

from ord_interface.editor.api.database import get_cursor

router = APIRouter(tags=["users"])


@router.get("/create_user")
async def create_user(user_name: str):
    """Creates a new user and returns the associated user ID."""
    user_id = uuid4().hex
    with get_cursor() as cursor:
        cursor.execute("INSERT INTO users (user_id, user_name) VALUES (%s, %s)", (user_id, user_name))
    return user_id


@router.get("/delete_user")
async def delete_user(user_id: str):
    """Deletes a user."""
    with get_cursor() as cursor:
        cursor.execute("DELETE FROM users WHERE user_id = %s", (user_id,))
