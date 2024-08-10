"""Celery configuration."""

import os

HOST = os.environ.get("REDIS_HOST", "localhost")
PORT = os.environ.get("REDIS_PORT", "6379")

# pylint: disable=invalid-name

broker_url = f"redis://{HOST}:{PORT}/1"
imports = ("search",)
result_backend = f"redis://{HOST}:{PORT}/2"
