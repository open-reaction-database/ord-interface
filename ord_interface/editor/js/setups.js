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

goog.module('ord.setups');
goog.module.declareLegacyNamespace();

const asserts = goog.require('goog.asserts');

const codes = goog.require('ord.codes');
const utils = goog.require('ord.utils');

const ReactionSetup = goog.require('proto.ord.ReactionSetup');
const ReactionEnvironment = goog.require('proto.ord.ReactionSetup.ReactionEnvironment');
const ReactionEnvironmentType = goog.require('proto.ord.ReactionSetup.ReactionEnvironment.ReactionEnvironmentType');
const Vessel = goog.require('proto.ord.Vessel');
const VesselAttachment = goog.require('proto.ord.VesselAttachment');
const VesselAttachmentType = goog.require('proto.ord.VesselAttachment.VesselAttachmentType');
const VesselMaterial = goog.require('proto.ord.VesselMaterial');
const VesselMaterialType = goog.require('proto.ord.VesselMaterial.VesselMaterialType');
const VesselPreparation = goog.require('proto.ord.VesselPreparation');
const VesselPreparationType = goog.require('proto.ord.VesselPreparation.VesselPreparationType');
const VesselType = goog.require('proto.ord.Vessel.VesselType');
const Volume = goog.require('proto.ord.Volume');

exports = {
  load,
  unload,
  addVesselPreparation,
  addVesselAttachment,
  validateSetup
};

/**
 * Adds and populates the reaction setup section in the form.
 * @param {!ReactionSetup} setup
 */
function load(setup) {
  loadVessel(setup.getVessel());
  const isAutomated = setup.hasIsAutomated() ? setup.getIsAutomated() : null;
  const setupAutomated = $('#setup_automated');
  utils.setOptionalBool(setupAutomated, isAutomated);
  if (isAutomated) {
    $('#automation_platform').show();
  }

  setupAutomated.on('change', function() {
    if (utils.getOptionalBool(setupAutomated)) {
      $('#automation_platform').show();
    } else {
      $('#automation_platform').hide();
    }
  });

  const platform = setup.getAutomationPlatform();
  $('#setup_platform').text(platform);

  const automationCodeMap = setup.getAutomationCodeMap();
  codes.load(automationCodeMap);

  const environment = setup.getEnvironment();
  if (environment != null) {
    utils.setSelector($('#setup_environment_type'), environment.getType());
    $('#setup_environment_details').text(environment.getDetails());
  }
}

/**
 * Adds and populates the reaction vessel section of the form.
 * @param {?Vessel} vessel
 */
function loadVessel(vessel) {
  if (!vessel) {
    return;
  }
  utils.setSelector($('#setup_vessel_type'), vessel.getType());
  $('#setup_vessel_details').text(vessel.getDetails());
  const material = vessel.getMaterial();
  if (material) {
    utils.setSelector($('#setup_vessel_material'), material.getType());
    $('#setup_vessel_material_details').text(material.getDetails());
  }
  const preparations = vessel.getPreparationsList();
  preparations.forEach(preparation => {
    const node = addVesselPreparation();
    utils.setSelector(
        $('.setup_vessel_preparation_type', node), preparation.getType());
    $('.setup_vessel_preparation_details', node).text(preparation.getDetails());
  });
  const attachments = vessel.getAttachmentsList();
  attachments.forEach(attachment => {
    const node = addVesselAttachment();
    utils.setSelector(
        $('.setup_vessel_attachment_type', node), attachment.getType());
    $('.setup_vessel_attachment_details', node).text(attachment.getDetails());
  });
  if (vessel.hasVolume()) {
    const volume = vessel.getVolume();
    utils.writeMetric('#setup_vessel_volume', volume);
  }
}

/**
 * Fetches the reaction setup from the form.
 * @return {!ReactionSetup}
 */
function unload() {
  const setup = new ReactionSetup();

  const vessel = unloadVessel();
  if (!utils.isEmptyMessage(vessel)) {
    setup.setVessel(vessel);
  }

  const isAutomated = utils.getOptionalBool($('#setup_automated'));
  if (isAutomated !== null) {
    setup.setIsAutomated(isAutomated);
  }

  setup.setAutomationPlatform(
      asserts.assertString($('#setup_platform').text()));

  const automationCodeMap = setup.getAutomationCodeMap();
  codes.unload(automationCodeMap);

  const environment = new ReactionEnvironment();
  const environmentType =
      utils.getSelectorText($('#setup_environment_type')[0]);
  environment.setType(ReactionEnvironmentType[environmentType]);
  environment.setDetails(
      asserts.assertString($('#setup_environment_details').text()));
  if (!utils.isEmptyMessage(environment)) {
    setup.setEnvironment(environment);
  }

  return setup;
}

/**
 * Fetches the reaction vessel information from the form.
 * @return {!Vessel}
 */
function unloadVessel() {
  const vessel = new Vessel();
  const vesselType = utils.getSelectorText($('#setup_vessel_type')[0]);
  vessel.setType(VesselType[vesselType]);
  vessel.setDetails(asserts.assertString($('#setup_vessel_details').text()));

  const material = new VesselMaterial();
  const materialType = utils.getSelectorText($('#setup_vessel_material')[0]);
  material.setType(VesselMaterialType[materialType]);
  material.setDetails(
      asserts.assertString($('#setup_vessel_material_details').text()));
  if (!utils.isEmptyMessage(material)) {
    vessel.setMaterial(material);
  }

  const preparations = [];
  $('.setup_vessel_preparation').each(function(index, node) {
    node = $(node);
    if (utils.isTemplateOrUndoBuffer(node)) {
      // The template.
      return;
    }
    const preparation = new VesselPreparation();
    const preparationType =
        utils.getSelectorText($('.setup_vessel_preparation_type', node)[0]);
    preparation.setType(VesselPreparationType[preparationType]);
    preparation.setDetails(asserts.assertString(
        $('.setup_vessel_preparation_details', node).text()));
    if (!utils.isEmptyMessage(preparation)) {
      preparations.push(preparation);
    }
  });
  vessel.setPreparationsList(preparations);

  const attachments = [];
  $('.setup_vessel_attachment').each(function(index, node) {
    node = $(node);
    if (utils.isTemplateOrUndoBuffer(node)) {
      return;
    }
    const attachment = new VesselAttachment();
    const attachmentType =
        utils.getSelectorText($('.setup_vessel_attachment_type', node)[0]);
    attachment.setType(VesselAttachmentType[attachmentType]);
    attachment.setDetails(asserts.assertString(
        $('.setup_vessel_attachment_details', node).text()));
    if (!utils.isEmptyMessage(attachment)) {
      attachments.push(attachment);
    }
  });
  vessel.setAttachmentsList(attachments);

  const volume = utils.readMetric('#setup_vessel_volume', new Volume());
  if (!utils.isEmptyMessage(volume)) {
    vessel.setVolume(volume);
  }

  return vessel;
}

/**
 * Adds a new vessel preparation section to the form.
 * @return {!jQuery} The node of the newly added div.
 */
function addVesselPreparation() {
  return utils.addSlowly(
      '#setup_vessel_preparation_template', $('#setup_vessel_preparations'));
}

/**
 * Adds a new vessel attachment section to the form.
 * @return {!jQuery} The node of the newly added div.
 */
function addVesselAttachment() {
  return utils.addSlowly(
      '#setup_vessel_attachment_template', $('#setup_vessel_attachments'));
}

/**
 * Validates the reaction setup defined in the form.
 * @param {!jQuery} node The node containing the reaction setup div.
 * @param {?jQuery=} validateNode The target div for validation results.
 */
function validateSetup(node, validateNode = null) {
  const setup = unload();
  utils.validate(setup, 'ReactionSetup', node, validateNode);
}
