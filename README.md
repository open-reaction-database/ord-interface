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

### `npm start`

Runs the app in the development mode.\
Open [http://localhost:3000](http://localhost:3000) to view it in your browser.

The page will reload when you make changes.\
You may also see any lint errors in the console.

### `npm test`

Launches the test runner in the interactive watch mode.\
See the section about [running tests](https://facebook.github.io/create-react-app/docs/running-tests) for more information.

### `npm run build`

Builds the app for production to the `build` folder.\
It correctly bundles React in production mode and optimizes the build for the best performance.

The build is minified and the filenames include the hashes.\
Your app is ready to be deployed!

See the section about [deployment](https://facebook.github.io/create-react-app/docs/deployment) for more information.