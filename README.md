# ord-interface

Web interface and FastAPI server for the [Open Reaction Database](https://open-reaction-database.org).

## Project layout

- **`app/`** — Vue single-page application (frontend)
  - `src/components/` — reusable Vue components
  - `src/views/` — routed pages (browse, search, dataset/reaction views)
  - `src/router/`, `src/styles/`, `src/utils/`, `src/assets/`
- **`ord_interface/`** — Python backend
  - `api/` — FastAPI server (search, view, queries)
  - `client/` — utilities for building the Postgres + rdkit-cartridge database
  - `visualization/` — Jinja filters and helpers for rendering reactions and molecules
  - `editor/` — legacy reaction-submission UI (being removed)

The Vue frontend talks to the FastAPI server, so both processes need to be running for the UI to work. On Windows, use WSL.

## Local development

Prerequisites: [`uv`](https://docs.astral.sh/uv/), Postgres with the [rdkit cartridge](https://www.rdkit.org/docs/Cartridge.html), [`node` + `npm`](https://nodejs.org/). Docker is optional, for the bundled full-stack image.

### 1. Install Python dependencies

```bash
git clone https://github.com/open-reaction-database/ord-interface
cd ord-interface
# Apple silicon: use `conda install postgresql` instead of the rdkit channel.
conda install -c rdkit rdkit-postgresql
uv sync
```

### 2. Run the API

#### Option A: Docker (bundled Postgres + nginx + gunicorn)

```shell
cd ord_interface
./build_test_database.sh
# Apple silicon: append `--build-arg="ARCH=aarch_64"`.
docker build --file Dockerfile -t openreactiondatabase/ord-interface ..
docker compose up
```

`docker compose up` binds Postgres on host port `5432` — stop any local Postgres first.

#### Option B: FastAPI dev server (in-process testing database)

```shell
cd ord_interface/api
ORD_INTERFACE_TESTING=TRUE fastapi dev main.py --port=5000
```

### 3. Run the Vue SPA

Download [Ketcher v2.5.1](https://github.com/epam/ketcher/releases/download/v2.5.1/ketcher-standalone-2.5.1.zip) and extract its contents into `./app/src/ketcher`, then:

```shell
cd app
npm install
npm run serve
```

Open <http://localhost:8080>. The page hot-reloads on edit.

## Deployment

Production deployment is managed via Pulumi/ECS in the [`ord-infrastructure`](https://github.com/open-reaction-database/ord-infrastructure) repository.

## Updating ord-schema

### Minor changes (e.g. enum additions)

- Update the `ord-schema` specifier in `pyproject.toml` under `[project] dependencies`, then run `uv lock` to refresh `uv.lock`.
- Update the matching version in [`ord_interface/Dockerfile`](./ord_interface/Dockerfile).
- Contact the ORD site administrator to roll out the new version to staging.

### Major changes (new message types, etc.)

The Vue app needs to be updated to surface new fields. Please discuss any such proposals with the ORD team before opening a PR.
