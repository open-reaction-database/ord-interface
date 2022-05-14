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

goog.module('ord.reaction');
goog.module.declareLegacyNamespace();

const asserts = goog.require('goog.asserts');

const conditions = goog.require('ord.conditions');
const electro = goog.require('ord.electro');
const flows = goog.require('ord.flows');
const identifiers = goog.require('ord.identifiers');
const illumination = goog.require('ord.illumination');
const inputs = goog.require('ord.inputs');
const notes = goog.require('ord.notes');
const observations = goog.require('ord.observations');
const outcomes = goog.require('ord.outcomes');
const pressure = goog.require('ord.pressure');
const provenance = goog.require('ord.provenance');
const setups = goog.require('ord.setups');
const stirring = goog.require('ord.stirring');
const temperature = goog.require('ord.temperature');
const uploads = goog.require('ord.uploads');
const utils = goog.require('ord.utils');
const workups = goog.require('ord.workups');

const Dataset = goog.require('proto.ord.Dataset');
const Reaction = goog.require('proto.ord.Reaction');

exports = {
  commit,
  downloadReaction,
  initFromDataset,
  initFromReactionId,
  unloadReaction,
  validateReaction,
};

const session = utils.session;

/**
 * Initializes the form.
 * @param {!Reaction} reaction Reaction proto to load.
 *
 * TODO(kearnes): Many undefined properties here.
 * @suppress {missingProperties}
 */
async function init(reaction) {
  // Initialize all the template popup menus.
  $('.selector').each((index, node) => utils.initSelector($(node)));
  $('.optional_bool').each((index, node) => utils.initOptionalBool($(node)));
  // Enable all the editable text fields.
  $('.edittext').attr('contentEditable', 'true');
  // Initialize all the validators.
  $('.validate').each((index, node) => utils.initValidateNode($(node)));
  // Initialize validation handlers that don't go in "add" methods.
  initValidateHandlers();
  // Initialize tooltips.
  $('[data-toggle="tooltip"]').tooltip();
  // Show "save" on modifications.
  utils.listen($('body'));
  // Load Ketcher content into an element with attribute role="application".
  document.getElementById('ketcher-iframe').contentWindow.ketcher.initKetcher();
  // Initialize the UI with the Reaction.
  loadReaction(reaction);
  let reactionId = reaction.getReactionId();
  if (!reactionId) {
    reactionId = 'Reaction ' + session.index;
  }
  $('#reaction_id').text(reactionId);
  utils.clean();
  // Trigger reaction-level validation.
  await validateReaction();
  // Initialize autosave being on.
  utils.toggleAutosave();
  // Signal to tests that the DOM is initialized.
  utils.ready();
}

/**
 * Initializes the form from a Dataset name and Reaction index.
 * @param {string} fileName Path to a Dataset proto.
 * @param {number} index The index of this Reaction in the Dataset.
 */
async function initFromDataset(fileName, index) {
  session.fileName = fileName;
  session.index = index;
  // Fetch the Dataset containing the Reaction proto.
  session.dataset = await utils.getDataset(fileName);
  asserts.assertInstanceof(session.dataset, Dataset);  // Type hint.
  const reaction = session.dataset.getReactionsList()[index];
  await init(reaction);
}

/**
 * Initializes the form from a Reaction ID.
 * @param {string} reactionId
 */
async function initFromReactionId(reactionId) {
  const reaction = await utils.getReactionById(reactionId);
  asserts.assertInstanceof(reaction, Reaction);  // Type hint.
  // NOTE(kearnes): Without this next line, `reaction` will be
  // partial/incomplete, and I have no idea why.
  console.log(reaction.toObject());
  await init(reaction);
  $('#dataset_context').hide();
}

/**
 * Updates the visual summary of the current reaction.
 * @param {!Reaction} reaction
 * @return {!Promise}
 */
function renderReaction(reaction) {
  return new Promise(resolve => {
    const xhr = new XMLHttpRequest();
    xhr.open('POST', '/render/reaction');
    const binary = reaction.serializeBinary();
    xhr.responseType = 'json';
    xhr.onload = function() {
      if (xhr.response !== null) {
        $('#reaction_render').html(asserts.assertString(xhr.response));
      }
      resolve();
    };
    xhr.send(binary);
  });
}

/**
 * Validates the current reaction.
 */
async function validateReaction() {
  const node = $('#sections');
  const validateNode = $('#reaction_validate');
  const reaction = unloadReaction();
  utils.validate(reaction, 'Reaction', node, validateNode);
  // Trigger all submessages to validate.
  $('.validate:visible:not(#reaction_validate)').trigger('click');
  // Render reaction as an HTML block.
  await renderReaction(reaction);
}

/**
 * Downloads the current reaction as a serialized Reaction proto.
 */
function downloadReaction() {
  const reaction = unloadReaction();
  const binary = reaction.serializeBinary();
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/reaction/download');
  xhr.onload = () => {
    // Make the browser write the file.
    const url = URL.createObjectURL(new Blob([xhr.response]));
    const link = document.createElement('a');
    link.setAttribute('href', url);
    link.setAttribute('download', 'reaction.pbtxt');
    document.body.appendChild(link);
    link.click();
  };
  xhr.send(binary);
}

/**
 * Adds and populates the form with the given reaction.
 * @param {!Reaction} reaction
 */
function loadReaction(reaction) {
  const identifiersList = reaction.getIdentifiersList();
  identifiers.load(identifiersList);
  const inputsMap = reaction.getInputsMap();
  // Reactions start with an input by default.
  if (inputsMap.getLength()) {
    inputs.load(inputsMap);
  } else {
    inputs.add($('#inputs'));
  }
  const setupMessage = reaction.getSetup();
  if (setupMessage) {
    setups.load(setupMessage);
  }
  const conditionsMessage = reaction.getConditions();
  if (conditionsMessage) {
    conditions.load(conditionsMessage);
  }
  const notesMessage = reaction.getNotes();
  if (notesMessage) {
    notes.load(notesMessage);
  }
  const observationsList = reaction.getObservationsList();
  observations.load(observationsList);

  const workupsList = reaction.getWorkupsList();
  workups.load(workupsList);

  const outcomesList = reaction.getOutcomesList();
  // Reactions start with an outcome by default.
  if (outcomesList.length) {
    outcomes.load(outcomesList);
  } else {
    outcomes.add();
  }

  const provenanceMessage = reaction.getProvenance();
  if (provenanceMessage) {
    provenance.load(provenanceMessage);
  }
  $('#reaction_id').text(reaction.getReactionId());

  // Clean up floating point entries.
  $('.floattext').each(function() {
    const node = $(this);
    if (node.text() !== '') {
      node.text(utils.prepareFloat(parseFloat(node.text())));
    }
  });
}

/**
 * Fetches the current reaction from the form.
 * @return {!Reaction}
 */
function unloadReaction() {
  const reaction = new Reaction();
  const identifiersList = identifiers.unload();
  reaction.setIdentifiersList(identifiersList);

  const inputsMap = reaction.getInputsMap();
  // isEmptyMessage check occurs in inputs.unload.
  inputs.unload(inputsMap);

  const setupMessage = setups.unload();
  if (!utils.isEmptyMessage(setupMessage)) {
    reaction.setSetup(setupMessage);
  }

  const conditionsMessage = conditions.unload();
  if (!utils.isEmptyMessage(conditionsMessage)) {
    reaction.setConditions(conditionsMessage);
  }

  const notesMessage = notes.unload();
  if (!utils.isEmptyMessage(notesMessage)) {
    reaction.setNotes(notesMessage);
  }

  const observationsList = observations.unload();
  reaction.setObservationsList(observationsList);

  const workupsList = workups.unload();
  reaction.setWorkupsList(workupsList);

  const outcomesList = outcomes.unload();
  reaction.setOutcomesList(outcomesList);

  const provenanceMessage = provenance.unload();
  if (!utils.isEmptyMessage(provenanceMessage)) {
    reaction.setProvenance(provenanceMessage);
  }

  // Setter does nothing when passed an empty string.
  reaction.setReactionId(asserts.assertString($('#reaction_id').text()));
  return reaction;
}

/**
 * Initializes the validation handlers. Some nodes are dynamically added or
 * removed; we add their validation handlers when the nodes themselves are
 * added. However, other nodes are always present in the HTML, and aren't
 * dynamically added nor removed. We add live validation to these nodes here.
 */
function initValidateHandlers() {
  // For setup
  const setupNode = $('#section_setup');
  utils.addChangeHandler(setupNode, () => {
    setups.validateSetup(setupNode);
  });

  // For conditions
  const conditionNode = $('#section_conditions');
  utils.addChangeHandler(conditionNode, () => {
    conditions.validateConditions(conditionNode);
  });

  // For temperature
  const temperatureNode = $('#section_conditions_temperature');
  utils.addChangeHandler(temperatureNode, () => {
    temperature.validateTemperature(temperatureNode);
  });

  // For pressure
  const pressureNode = $('#section_conditions_pressure');
  utils.addChangeHandler(pressureNode, () => {
    pressure.validatePressure(pressureNode);
  });

  // For stirring
  const stirringNode = $('#section_conditions_stirring');
  utils.addChangeHandler(stirringNode, () => {
    stirring.validateStirring(stirringNode);
  });

  // For illumination
  const illuminationNode = $('#section_conditions_illumination');
  utils.addChangeHandler(illuminationNode, () => {
    illumination.validateIllumination(illuminationNode);
  });

  // For electro
  const electroNode = $('#section_conditions_electro');
  utils.addChangeHandler(electroNode, () => {
    electro.validateElectro(electroNode);
  });

  // For flow
  const flowNode = $('#section_conditions_flow');
  utils.addChangeHandler(flowNode, () => {
    flows.validateFlow(flowNode);
  });

  // For notes
  const notesNode = $('#section_notes');
  utils.addChangeHandler(notesNode, () => {
    notes.validateNotes(notesNode);
  });

  // For provenance
  const provenanceNode = $('#section_provenance');
  utils.addChangeHandler(provenanceNode, () => {
    provenance.validateProvenance(provenanceNode);
  });
}

/**
 * Writes the current reaction to disk.
 */
function commit() {
  if (!session.dataset) {
    // Do nothing when there is no Dataset; e.g. when viewing reactions by ID.
    return;
  }
  const reaction = unloadReaction();
  asserts.assertInstanceof(session.dataset, Dataset);  // Type hint.
  const reactions = session.dataset.getReactionsList();
  reactions[asserts.assertNumber(session.index)] = reaction;
  const fileName = asserts.assertString(session.fileName);
  utils.putDataset(fileName, session.dataset);
  uploads.putAll(fileName);
}
