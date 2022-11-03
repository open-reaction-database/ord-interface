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

## React App

```shell
cd app
```
In the `./app` directory, you can run:

### `npm run serve`

Runs the app in the development mode.\
Open [http://localhost:8080](http://localhost:8080) to view it in your browser.

The page will reload when you make changes.\
You may also see any lint errors in the console.

### `npm run build`

Builds the app for production to the `build` folder.\

The build is minified and the filenames include the hashes.\
Your app is ready to be deployed!

See the section about [deployment](https://facebook.github.io/create-react-app/docs/deployment) for more information.