# ord-interface

Backend and frontend code for the [ORD search interface](https://client.open-reaction-database.org/).

## Installation

```shell
$ git clone git@github.com:open-reaction-database/ord-interface.git
$ cd ord-interface
$ pip install -e .
```

To build and launch the interface (available at `http://localhost:5001`):

```shell
$ cd ord_interface
$ ./build_test_database.sh
# If you are running on Apple silicon, append `--build-arg="ARCH=aarch_64"` to the next command.
$ docker build --file Dockerfile -t openreactiondatabase/ord-interface ..
$ docker compose up
```

## Development

To start a Flask server in development mode:

```shell
$ cd ord_interface
$ ./build_test_database.sh
# Start the database backend.
$ docker run -d -p 5432:5432 openreactiondatabase/ord-postgres:test
# Start the development server.
$ POSTGRES_USER=postgres POSTGRES_PASSWORD=postgres FLASK_APP=interface.py FLASK_ENV=development python -m flask run
```
