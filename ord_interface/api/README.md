# Open Reaction Database API

## Local development

Set up your environment by installing postgres with the RDKit cartridge:

```shell
conda install postgresql python=3.10
# Note that rdkit-postgresql is currently only available on x86_64.
conda install -c rdkit rdkit-postgresql
cd ord-interface
pip install -e ".[tests]"
```

To start a test server locally, run:

```shell
cd ord-interface/ord_interface/api
ORD_INTERFACE_TESTING=TRUE fastapi dev main.py
```

This does a few things:
* Creates a test postgres database (using `testing.postgresql`)
* Populates the database using the datasets in `testdata`
* Launches a FastAPI server (see https://fastapi.tiangolo.com/#run-it)

To view the API docs, navigate to http://localhost:8000/docs.

To test with the Vue interface:

```shell
# In one terminal session:
cd ord-interface/ord_interface/api
ORD_INTERFACE_TESTING=TRUE fastapi dev main --port=5000
# In a second terminal session:
cd ord-interface/app
npm run serve
```