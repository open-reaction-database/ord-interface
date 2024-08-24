# ORD Editor API

## Local development

1. Install dependencies

   * PostgreSQL
   * Python >= 3.10

2. Install the package

    ```shell
    git clone git@github.com:open-reaction-database/ord-interface.git
    cd ord-interface
    pip install -e ".[tests]"
    ```

3. Run the FastAPI server

    ```shell
    cd ord_interface/editor/api
    ORD_EDITOR_TESTING=TRUE fastapi dev main.py
    ```
    
    This creates a test PostgreSQL database and starts the server at http://localhost:8000. Navigate to
    http://localhost:8000/docs for the interactive Swagger docs.
