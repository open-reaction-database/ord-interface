#!/bin/bash
#
# Runs the app with gunicorn behind a nginx proxy.
set -e

# Start nginx server.
nginx -g 'daemon off;' &

# Start gunicorn.
LOG_FORMAT='GUNICORN %(t)s %({user-id}o)s %(U)s %(s)s %(L)s %(b)s %(f)s "%(r)s" "%(a)s"'
gunicorn ord_interface.interface:app \
  --bind unix:/run/gunicorn.sock \
  --workers 2 \
  --access-logfile - \
  --access-logformat "${LOG_FORMAT}"
