# Open Reaction Database API

## Local development

Set up your environment by installing postgres with the RDKit cartridge:

```shell
conda install -c rdkit rdkit-postgresql
```

To start a test server locally, run:

```shell
ORD_INTERFACE_TESTING=TRUE fastapi dev main.py
```

This does a few things:
* Creates a test postgres database (using `testing.postgresql`)
* Populates the database using the datasets in `testdata`
* Launches a FastAPI server (see https://fastapi.tiangolo.com/#run-it)

To view the API docs, navigate to http://localhost:8000/docs.
