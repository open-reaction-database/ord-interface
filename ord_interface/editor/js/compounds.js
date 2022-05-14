/**
 * Copyright 2020 Open Reaction Database Project Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

goog.module('ord.compounds');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  unloadCompound,
  loadIntoCompound,
  add,
  validateCompound,
  drawIdentifier,
  addFeature,
  addNameIdentifier,
  addIdentifier,
  addPreparation,
  loadIdentifier,
  loadFeature,
  unloadFeature,
  unloadIdentifiers,
  renderCompound,
};

const asserts = goog.require('goog.asserts');

const JspbMap = goog.requireType('jspb.Map');

const amounts = goog.require('ord.amounts');
const data = goog.require('ord.data');
const uploads = goog.require('ord.uploads');
const utils = goog.require('ord.utils');

const Compound = goog.require('proto.ord.Compound');
const CompoundIdentifier = goog.require('proto.ord.CompoundIdentifier');
const IdentifierType = goog.require('proto.ord.CompoundIdentifier.IdentifierType');
const CompoundPreparation = goog.require('proto.ord.CompoundPreparation');
const PreparationType = goog.require('proto.ord.CompoundPreparation.PreparationType');
const Source = goog.require('proto.ord.Compound.Source');
const Data = goog.require('proto.ord.Data');
const ProductCompound = goog.require('proto.ord.ProductCompound');
const ReactionRoleType = goog.require('proto.ord.ReactionRole.ReactionRoleType');

/**
 * Adds and populates the form's fields describing multiple compounds for a
 * single reaction input.
 * @param {!jQuery} node The div corresponding to the reaction input to which
 *     compound definitions should be added.
 * @param {!Array<!Compound>} compounds
 */
function load(node, compounds) {
  compounds.forEach(compound => loadCompound(node, compound));
}

/**
 * Adds fields describing a new component to an existing reaction input in the
 * form and populates them according to the provided compound.
 * @param {!jQuery} root The div corresponding to the reaction input to which
 *     a new compound definition should be added.
 * @param {!Compound} compound
 */
function loadCompound(root, compound) {
  const node = add(root);
  loadIntoCompound(node, compound);
}

/**
 * Adds and populates the form's fields describing a compound.
 * @param {!jQuery} node The div corresponding to the compound whose fields
 *     should be updated.
 * @param {?Compound} compound
 */
function loadIntoCompound(node, compound) {
  if (!compound) {
    return;
  }
  const reactionRole = compound.getReactionRole();
  utils.setSelector($('.component_reaction_role', node), reactionRole);
  $('.component_reaction_role', node).trigger('change');

  const isLimiting = compound.hasIsLimiting() ? compound.getIsLimiting() : null;
  utils.setOptionalBool($('.component_limiting', node), isLimiting);

  const identifiers = compound.getIdentifiersList();
  identifiers.forEach(identifier => loadIdentifier(node, identifier));

  const amount = compound.getAmount();
  amounts.load(node, amount);

  const preparations = compound.getPreparationsList();
  preparations.forEach(preparation => {
    const preparationNode = addPreparation(node);
    loadPreparation(preparationNode, preparation);
  });
  if (compound.hasSource()) {
    const source = compound.getSource();
    loadSource(node, source);
  }
  const features = compound.getFeaturesMap();
  features.forEach(function(feature, name) {
    const featureNode = addFeature(node);
    loadFeature(featureNode, name, feature);
  });
}

/**
 * Adds fields describing a new identifier to an existing compound in the form
 * and populates them according to the provided identifier.
 * @param {!jQuery} compoundNode The div corresponding to the compound to which
 *     a new compound definition should be added.
 * @param {!CompoundIdentifier} identifier
 */
function loadIdentifier(compoundNode, identifier) {
  const node = addIdentifier(compoundNode);
  const value = identifier.getValue();
  $('.component_identifier_value', node).first().text(value);
  utils.setSelector(node, identifier.getType());
  $('.component_identifier_details', node)
      .first()
      .text(identifier.getDetails());
}

/**
 * Adds and populates the form's fields describing a compound preparation.
 * @param {!jQuery} node The div corresponding to the preparation that should be
 *     updated on the form.
 * @param {!CompoundPreparation} preparation
 */
function loadPreparation(node, preparation) {
  const type = preparation.getType();
  utils.setSelector($('.component_compound_preparation_type', node), type);
  const details = preparation.getDetails();
  $('.component_compound_preparation_details', node).text(details);
  const reaction = preparation.getReactionId();
  $('.component_compound_preparation_reaction', node).text(reaction);
}

/**
 * Adds and populates the form's fields describing a compound's source.
 * @param {!jQuery} compoundNode The div corresponding to the compound whose
 *     source information should be updated on the form.
 * @param {?Source} source
 */
function loadSource(compoundNode, source) {
  if (!source) {
    return;
  }
  const node = $('fieldset.source', compoundNode);
  $('.component_source_vendor', node).text(source.getVendor());
  $('.component_source_id', node).text(source.getId());
  $('.component_source_lot', node).text(source.getLot());
}

/**
 * Reads and returns a list of compounds defined within part of the form.
 * @param {!jQuery} node The div corresponding to the reaction inputs whose
 *     compounds should be read from the form.
 * @return {!Array<!Compound>}
 */
function unload(node) {
  const compounds = [];
  $('.component', node).each(function(index, compoundNode) {
    compoundNode = $(compoundNode);
    if (!compoundNode.attr('id')) {
      // Not a template.
      const compound = unloadCompound(compoundNode);
      if (!utils.isEmptyMessage(compound)) {
        compounds.push(compound);
      }
    }
  });
  return compounds;
}

/**
 * Reads and returns a single compound as defined on the form.
 * @param {!jQuery} node The div corresponding to the compound whose definition
 *     should be read from the form.
 * @return {!Compound}
 */
function unloadCompound(node) {
  const compound = new Compound();

  const reactionRole =
      utils.getSelectorText($('.component_reaction_role', node)[0]);
  compound.setReactionRole(ReactionRoleType[reactionRole]);

  // Only call setIsLimiting if this is a reactant Compound.
  if (utils.getSelectorText($('.component_reaction_role', node)[0]) ===
      'REACTANT') {
    const isLimiting = utils.getOptionalBool($('.component_limiting', node));
    if (isLimiting !== null) {
      compound.setIsLimiting(isLimiting);
    }
  }

  const identifiers = unloadIdentifiers(node);
  if (identifiers.some(e => !utils.isEmptyMessage(e))) {
    compound.setIdentifiersList(identifiers);
  }

  const amount = amounts.unload(node);
  if (!utils.isEmptyMessage(amount)) {
    compound.setAmount(amount);
  }

  const preparations = [];
  $('.component_preparation', node).each(function(index, preparationNode) {
    asserts.assertInstanceof(preparationNode, Element);
    const preparation = unloadPreparation(preparationNode);
    preparations.push(preparation);
  });
  if (preparations.some(e => !utils.isEmptyMessage(e))) {
    compound.setPreparationsList(preparations);
  }

  const source = unloadSource(node);
  if (!utils.isEmptyMessage(source)) {
    compound.setSource(source);
  }

  const featuresMap = compound.getFeaturesMap();
  $('.feature', node).each(function(index, featureNode) {
    featureNode = $(featureNode);
    if (!featureNode.attr('id')) {
      unloadFeature(featureNode, featuresMap);
    }
  });

  return compound;
}

/**
 * Reads and returns a list of compound identifiers for a single compound as
 * defined on the form.
 * @param {!jQuery} node The div corresponding to the compound whose identifiers
 *     should be read from the form.
 * @return {!Array<!CompoundIdentifier>}
 */
function unloadIdentifiers(node) {
  const identifiers = [];
  $('.component_identifier', node).each(function(index, node) {
    node = $(node);
    if (!utils.isTemplateOrUndoBuffer(node)) {
      const identifier = unloadIdentifier(node);
      if (!utils.isEmptyMessage(identifier)) {
        identifiers.push(identifier);
      }
    }
  });
  return identifiers;
}

/**
 * Reads and returns a single compound identifier as defined on the form.
 * @param {!jQuery} node The div corresponding to the compound identifier that
 *     should be read from the form.
 * @return {!CompoundIdentifier}
 */
function unloadIdentifier(node) {
  const identifier = new CompoundIdentifier();

  const value = $('.component_identifier_value', node).text();
  identifier.setValue(asserts.assertString(value));
  const type = utils.getSelectorText(node[0]);
  identifier.setType(IdentifierType[type]);
  const details = $('.component_identifier_details', node).text();
  identifier.setDetails(asserts.assertString(details));
  return identifier;
}

/**
 * Reads and returns a single compound preparation as defined on the form.
 * @param {!Element} node The div corresponding to a compound preparation that
 *     should be read from the form.
 * @return {!CompoundPreparation}
 */
function unloadPreparation(node) {
  const preparation = new CompoundPreparation();
  const type =
      utils.getSelectorText($('.component_compound_preparation_type', node)[0]);
  preparation.setType(PreparationType[type]);
  const details = $('.component_compound_preparation_details', node).text();
  preparation.setDetails(asserts.assertString(details));
  const reaction = $('.component_compound_preparation_reaction', node).text();
  preparation.setReactionId(asserts.assertString(reaction));
  return preparation;
}

/**
 * Sets the source information fields of a compound according to the form.
 * @param {!jQuery} node The div corresponding to the compound whose source
 *     information should be read from the form.
 * @return {!Source}
 */
function unloadSource(node) {
  const source = new Source();
  const vendor = $('.component_source_vendor', node).text();
  source.setVendor(asserts.assertString(vendor));
  const lot = $('.component_source_lot', node).text();
  source.setLot(asserts.assertString(lot));
  const id = $('.component_source_id', node).text();
  source.setId(asserts.assertString(id));
  return source;
}

/**
 * Adds fields to the form corresponding to a new, empty compound definition as
 * specified by the component template with ID "component_template".
 * @param {!jQuery} root The div within which the new compound should be added.
 * @return {!jQuery} The node of the new component div.
 */
function add(root) {
  const node = utils.addSlowly('#component_template', $('.components', root));

  // Connect reaction role selection to limiting reactant field.
  const roleSelector = $('.component_reaction_role', node);
  roleSelector.on('change', function() {
    if (utils.getSelectorText(this) === 'REACTANT') {
      $('.limiting_reactant', node).show();
    } else {
      $('.limiting_reactant', node).hide();
    }
  });

  amounts.init(node);

  // Add live validation handling.
  utils.addChangeHandler(node, () => {
    validateCompound(node);
  });

  return node;
}

/**
 * Adds fields to the form corresponding to a new, empty compound identifier as
 * specified by the component identifier template with ID
 * "component_identifier_template".
 * @param {!jQuery} node The div within which the new identifier should be
 *     added.
 * @return {!jQuery} The node of the new compound identifier div.
 */
function addIdentifier(node) {
  const identifierNode = utils.addSlowly(
      '#component_identifier_template', $('.identifiers', node).first());

  const uploadButton = $('.component_identifier_upload', identifierNode);
  uploadButton.on('change', function() {
    if ($(this).is(':checked')) {
      $('.uploader', identifierNode).show();
      $('.component_identifier_value', identifierNode).hide();
      $('.text_upload', identifierNode).hide();
    } else {
      $('.uploader', identifierNode).hide();
      $('.component_identifier_value', identifierNode).show();
    }
  });
  uploads.initialize(identifierNode);
  return identifierNode;
}

/**
 * Adds new compound identifier(s) after prompting the user for the name of
 * a compound to add. A NAME-type identifier is always included. A SMILES-type
 * identifier is included when the name can be parsed.
 * @param {!jQuery} node The div corresponding to the compound to which the new
 *     identifiers should be added.
 */
function addNameIdentifier(node) {
  const name = prompt('Compound name: ');
  if (!(name)) {
    return;
  }
  const identifier = new CompoundIdentifier();
  identifier.setValue(name);
  identifier.setType(IdentifierType.NAME);
  loadIdentifier(node, identifier);

  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/resolve/name');
  xhr.responseType = 'json';
  xhr.onload = function() {
    if (xhr.response) {
      const smiles = xhr.response[0];
      const resolver = xhr.response[1];
      const identifier = new CompoundIdentifier();
      identifier.setValue(smiles);
      identifier.setType(IdentifierType.SMILES);
      identifier.setDetails('NAME resolved by the ' + resolver);
      loadIdentifier(node, identifier);
    }
    validateCompound(node);
  };
  xhr.send(name);
}

/**
 * Displays the Ketcher drawing tool for defining compound identifiers. A
 * callback is defined to create SMILES and MOLBLOCK-type identifiers; this
 * callback is triggered upon submission of a molecular drawing from the
 * Ketcher window.
 * @param {!jQuery} node The div corresponding to the compound to which the new
 *     identifiers should be added.
 *
 * NOTE(kearnes): Lots of undefined properties on `ketcher`.
 * NOTE(kearnes): `.modal()` is defined by jQuery Modal.
 * @suppress {missingProperties}
 */
function drawIdentifier(node) {
  // Get a reference to Ketcher, and to look nice, clear any old drawings.
  const ketcher =
      document.getElementById('ketcher-iframe').contentWindow.ketcher;
  ketcher.editor.struct(null);
  // Start the loading spinner.
  $('#ketcher-spinner').fadeIn(0);

  // First, pack the current Compound into a message.
  const compound = new Compound();
  const identifiers = unloadIdentifiers(node);
  if (identifiers.some(e => !utils.isEmptyMessage(e))) {
    compound.setIdentifiersList(identifiers);
  }
  // Then, try to resolve compound into a MolBlock.
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/ketcher/molfile');
  const binary = compound.serializeBinary();
  xhr.responseType = 'json';
  xhr.onload = function() {
    if (xhr.status === 200) {
      const molblock = xhr.response;
      // Set the molecule in ketcher.
      // Note: In case async / callback issues prove difficult,
      // a cleaner fix may be to put this entire xhr in a modal callback, then
      // toggle the modal.

      // If the modal is already open, we can simply set the molecule.
      const ketcherModal = $('#ketcher_modal');
      if (ketcherModal.hasClass('show')) {
        ketcher.setMolecule(molblock);
      }
      // Otherwise, we need to set up a callback, so that the molecule is set
      // only when Ketcher is open. (to prevent a graphical glitch)
      else {
        ketcherModal.on('shown.bs.modal', function() {
          // This callback should only be ever run once, so make sure to remove
          // it.
          ketcherModal.off('shown.bs.modal');
          ketcher.setMolecule(molblock);
        });
      }
    }
    // Now that we're done with (trying to) loading the molecule, hide the
    // spinner.
    $('#ketcher-spinner').fadeOut();
  };
  xhr.send(binary);
  // Finally, open the ketcher modal.
  $('#ketcher_modal').modal('show');
  // Define a callback so that when a user is done drawing, the new SMILES
  // string gets saved.
  ketcher.successCallback = function() {
    // Check if an existing SMILES/MolBlock identifier exists. If yes, remove.
    $('.component_identifier', node).each(function(index, node) {
      node = $(node);
      if (!utils.isTemplateOrUndoBuffer(node)) {
        const identifier = unloadIdentifier(node);
        if ((identifier.getType() === IdentifierType.SMILES) ||
            (identifier.getType() === IdentifierType.MOLBLOCK)) {
          utils.removeSlowly(node, '.component_identifier');
        }
      }
    });
    // Create new identifiers.
    if (ketcher.getSmiles()) {
      const xhr = new XMLHttpRequest();
      xhr.open('POST', '/canonicalize');
      xhr.responseType = 'json';
      xhr.onload = function() {
        const smilesIdentifier = new CompoundIdentifier();
        smilesIdentifier.setType(IdentifierType.SMILES);
        smilesIdentifier.setValue(asserts.assertString(xhr.response));
        smilesIdentifier.setDetails('Drawn with Ketcher');
        loadIdentifier(node, smilesIdentifier);
      };
      xhr.send(ketcher.getSmiles());
      const molfileIdentifier = new CompoundIdentifier();
      molfileIdentifier.setType(IdentifierType.MOLBLOCK);
      molfileIdentifier.setValue(ketcher.getMolfile());
      molfileIdentifier.setDetails('Drawn with Ketcher');
      loadIdentifier(node, molfileIdentifier);
    }
    validateCompound(node);
  };
}

/**
 * Adds fields to the form corresponding to a new, empty compound preparation
 * as specified by the component preparation template with ID
 * "component_preparation_template".
 * @param {!jQuery} node The div corresponding to the compound to which the new
 *     preparation should be added.
 * @return {!jQuery} The div corresponding to the new compound preparation.
 */
function addPreparation(node) {
  const PreparationNode = utils.addSlowly(
      '#component_preparation_template', $('.preparations', node));

  const typeSelector =
      $('.component_compound_preparation_type', PreparationNode);
  typeSelector.on('change', function() {
    if (utils.getSelectorText(this) === 'SYNTHESIZED') {
      $('.component_compound_preparation_reaction_id', PreparationNode)
          .css('display', 'inline-block');
    } else {
      $('.component_compound_preparation_reaction_id', PreparationNode)
          .css('display', 'none');
    }
  });

  return PreparationNode;
}

/**
 * Updates a png rendering of a compound as defined by its identifiers.
 * @param {!jQuery} node The div corresponding to the compound whose rendering
 *     should be updated.
 * @param {!Compound|!ProductCompound} compound
 */
function renderCompound(node, compound) {
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/render/compound');
  const binary = compound.serializeBinary();
  xhr.responseType = 'json';
  xhr.onload = function() {
    if (xhr.response) {
      $('.component_rendering', node).html(asserts.assertString(xhr.response));
    } else {
      $('.component_rendering', node).html('');
    }
  };
  xhr.send(binary);
}

/**
 * Validates the definition of a compound and updates the validation error
 * display node.
 * @param {!jQuery} node The div corresponding to the compound that should be
 *     read from the form and validated.
 * @param {?jQuery=} validateNode The div that is used to show the results of
 *     validation (i.e., success or errors).
 */
function validateCompound(node, validateNode = null) {
  const compound = unloadCompound(node);
  utils.validate(compound, 'Compound', node, validateNode);

  // Try to resolve compound structural identifiers. This is tied to
  // validation so the same trigger is used and we only have to unload the
  // compound once per update.
  renderCompound(node, compound);
}

/**
 * Adds a new feature section to the form.
 * @param {!jQuery} node Parent component node.
 * @return {!jQuery} The newly added parent node for the Data record.
 */
function addFeature(node) {
  const featureNode =
      utils.addSlowly('#feature_template', $('.features', node));
  data.addData(featureNode);
  return featureNode;
}

/**
 * Adds and populates a feature section in a Compound.
 * @param {!jQuery} node Parent component node.
 * @param {string} name The name of this Data record.
 * @param {!Data} feature
 */
function loadFeature(node, name, feature) {
  $('.feature_name', node).text(name);
  data.loadData(node, feature);
}

/**
 * Fetches a feature record defined in the form and adds it to `featuresMap`.
 * @param {!jQuery} node Root node for the Data record.
 * @param {!JspbMap<string, !Data>} featuresMap
 */
function unloadFeature(node, featuresMap) {
  const name = $('.feature_name', node).text();
  const dataMessage = data.unloadData(node);
  if (name || !utils.isEmptyMessage(dataMessage)) {
    featuresMap.set(asserts.assertString(name), dataMessage);
  }
}
