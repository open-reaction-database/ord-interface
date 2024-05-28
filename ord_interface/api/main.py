"""Open Reaction Database API."""

from fastapi import FastAPI

from ord_interface.api import client, editor

app = FastAPI()
app.include_router(client.router)
app.include_router(editor.router)


@app.get("/ketcher")
def ketcher():
    return app.redirect("/standalone/index.html")
