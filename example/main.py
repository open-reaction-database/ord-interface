"""Minimal flask example for ORD protocol buffers."""
from flask import Flask, jsonify, make_response, render_template, request
from ord_schema.proto import reaction_pb2
from ord_schema.validations import ValidationOptions, validate_message

app = Flask(__name__, template_folder=".")


@app.route("/")
def index():
    return render_template("main.html")


@app.route("/fetchReaction")
def fetchReaction():
    reaction = reaction_pb2.Reaction()
    reaction.identifiers.add(value="C(C)Cl.Br>>C(C)Br.Cl", type="REACTION_SMILES")
    response = make_response(reaction.SerializeToString())
    response.headers.set("Content-Type", "application/protobuf")
    return response


@app.route("/validateReaction", methods=["POST"])
def validate_reaction():
    reaction = reaction_pb2.Reaction()
    reaction.ParseFromString(request.get_data())
    options = ValidationOptions(require_provenance=True)
    output = validate_message(reaction, raise_on_error=False, options=options)
    return jsonify({"errors": output.errors, "warnings": output.warnings})
