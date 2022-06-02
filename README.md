# ord-interface

Backend and frontend code for the [ORD search interface](https://client.open-reaction-database.org/).

# Development setup

```bash
conda create -n ord-interface python==3.9
conda activate ord-interface
pip install -r requirements.txt
pip install .

cd ord_interface
./build_test_database.sh
docker build --file Dockerfile -t openreactiondatabase/ord-interface .. && docker compose up


```