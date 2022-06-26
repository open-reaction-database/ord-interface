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
$ docker build --file Dockerfile -t openreactiondatabase/ord-interface ..
$ docker compose up
```
