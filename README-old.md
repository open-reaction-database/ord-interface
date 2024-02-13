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

## Vue App

```shell
cd app
```
In the `./app` directory, you can run:

### Initialize app `npm i`

Install the necessary npm packages for the vue app.
Download [https://github.com/epam/ketcher/releases/download/v2.5.1/ketcher-standalone-2.5.1.zip](Ketcher) and extract the files into `./app/src/ketcher/`

### `npm run serve`

Runs the app in the development mode.
Open [http://localhost:8080](http://localhost:8080) to view it in your browser.
The vue app depends on the flask app api running on port 5000 (see the Development instructions for the flask app above)

The page will reload when you make changes.
You may also see any lint errors in the console.

### `npm run build`

Builds the app for production to the `dist` folder.\