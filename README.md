# ord-interface

Backend and frontend code for the [ORD search interface](https://client.open-reaction-database.org/).

## Installation

```shell
git clone git@github.com:open-reaction-database/ord-interface.git
cd ord-interface
pip install -e .
```

To build and launch the interface (available at `http://localhost:5001`):

```shell
cd ord_interface
./build_test_database.sh
docker build --file Dockerfile -t openreactiondatabase/ord-interface ..
docker compose up
```

## Development

### Search/browse interface

To start a Flask server in development mode:

```shell
cd ord_interface
./build_test_database.sh
# Start the database backend.
docker run -d -p 5432:5432 openreactiondatabase/ord-postgres:test
# Start the development server.
POSTGRES_USER=postgres POSTGRES_PASSWORD=postgres FLASK_APP=interface.py FLASK_ENV=development python -m flask run
```

### Editor

Because the closure libraries need to be rebuilt after every change, it is recommended to use the normal docker compose
path:

```shell
cd ord_interface
./build_test_database.sh  # This only needs to be run once (unless the test data are changed).
docker build --file Dockerfile -t openreactiondatabase/ord-interface ..
docker compose up
```
