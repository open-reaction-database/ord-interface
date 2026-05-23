# ord-interface

Web interface and FastAPI server for the [Open Reaction Database](https://open-reaction-database.org).

## Project layout

- **`app/`** — React single-page application (frontend), built with Vite
  - `src/components/` — reusable React components
  - `src/views/` — routed pages (browse, search, dataset/reaction views)
  - `src/router/`, `src/styles/`, `src/utils/`, `src/assets/`
- **`ord_interface/`** — Python backend
  - `api/` — FastAPI server (search, view, queries)
  - `client/` — utilities for building the Postgres + rdkit-cartridge database
  - `visualization/` — Jinja filters and helpers for rendering reactions and molecules
  - `editor/` — legacy reaction-submission UI (being removed)

The React frontend talks to the FastAPI server, so both processes need to be running for the UI to work. On Windows, use WSL.

## Local development

Prerequisites:

- [`uv`](https://docs.astral.sh/uv/)
- Postgres with the [rdkit cartridge](https://www.rdkit.org/docs/Cartridge.html)
- [`node` + `npm`](https://nodejs.org/)
- [`redis`](https://redis.io/) — the API uses it for search-task state
- Docker (optional, for the bundled full-stack image)

### 1. Install Python dependencies

```bash
git clone https://github.com/open-reaction-database/ord-interface
cd ord-interface
# Apple silicon: use `conda install postgresql` instead of the rdkit channel.
conda install -c rdkit rdkit-postgresql
uv sync
```

The `Local` option below shells out to `initdb`/`postgres` to bring up an
in-process test database, so the conda environment that has
`rdkit-postgresql` needs to stay activated (or its `bin/` on `PATH`) for the
rest of the steps.

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

#### Option B: Local (in-process testing database)

In your local conda environment
```shell
# Redis (the API uses it for search-task state); leaves it daemonized.
redis-server --daemonize yes

# Backend. ORD_INTERFACE_TESTING=TRUE spins up an in-process Postgres
# via testing.postgresql and loads the bundled test datasets.
ORD_INTERFACE_TESTING=TRUE uv run uvicorn ord_interface.api.main:app --port 5000 --reload
```

> **macOS note:** Apple's AirPlay Receiver service also listens on `*:5000`. If
> requests hit the wrong listener you'll get `HTTP 403` with a `Server:
> AirTunes` header. Either disable AirPlay Receiver (System Settings → General
> → AirDrop & Handoff), or bind uvicorn to a different port and update the
> proxy target in `app/vite.config.js`. Uvicorn binds `127.0.0.1` only, so
> `http://127.0.0.1:5000/...` will reach uvicorn regardless.

### 3. Run the React SPA

Download [Ketcher v2.5.1](https://github.com/epam/ketcher/releases/download/v2.5.1/ketcher-standalone-2.5.1.zip), extract it, and place the `standalone/` directory's contents inside `./app/public/ketcher/` (so that `app/public/ketcher/static/` and `app/public/ketcher/asset-manifest.json` exist). Then:

```shell
cd app
npm install
npm run serve
```

Open <http://localhost:8080>. Vite hot-reloads source changes; `npm run dev` and `npm run preview` are also available as aliases for `vite` and `vite preview` respectively.

## Deployment

Production deployment is managed via Pulumi/ECS in the [`ord-infrastructure`](https://github.com/open-reaction-database/ord-infrastructure) repository.

## Updating ord-schema

### Minor changes (e.g. enum additions)

- Update the `ord-schema` specifier in `pyproject.toml` under `[project] dependencies`, then run `uv lock` to refresh `uv.lock`.
- Update the matching version in [`ord_interface/Dockerfile`](./ord_interface/Dockerfile).
- Contact the ORD site administrator to roll out the new version to staging.

### Major changes (new message types, etc.)

The React app needs to be updated to surface new fields. Please discuss any such proposals with the ORD team before opening a PR.
