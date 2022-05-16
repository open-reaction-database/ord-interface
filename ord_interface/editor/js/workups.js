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

goog.module('ord.workups');
goog.module.declareLegacyNamespace();

const asserts = goog.require('goog.asserts');

const amounts = goog.require('ord.amounts');
const inputs = goog.require('ord.inputs');
const stirring = goog.require('ord.stirring');
const temperature = goog.require('ord.temperature');
const utils = goog.require('ord.utils');

const ReactionWorkup = goog.require('proto.ord.ReactionWorkup');
const WorkupType = goog.require('proto.ord.ReactionWorkup.WorkupType');
const Temperature = goog.require('proto.ord.Temperature');
const TemperatureConditions = goog.require('proto.ord.TemperatureConditions');
const Measurement = goog.require('proto.ord.TemperatureConditions.Measurement');
const MeasurementType = goog.require('proto.ord.TemperatureConditions.Measurement.MeasurementType');
const Time = goog.require('proto.ord.Time');

exports = {
  load,
  unload,
  add,
  addMeasurement,
  validateWorkup
};

/**
 * Adds and populates the reaction workup sections in the form.
 * @param {!Array<!ReactionWorkup>} workups
 */
function load(workups) {
  workups.forEach(workup => loadWorkup(workup));
}

/**
 * Adds and populates a reaction workup section in the form.
 * @param {!ReactionWorkup} workup
 */
function loadWorkup(workup) {
  const node = add();
  utils.setSelector($('.workup_type', node), workup.getType());
  $('.workup_type select', node).trigger('change');
  $('.workup_details', node).text(workup.getDetails());
  const duration = workup.getDuration();
  if (duration) {
    utils.writeMetric('.workup_duration', duration, node);
  }

  const input = workup.getInput();
  inputs.loadInputUnnamed($('.workup_input', node), input);

  const amount = workup.getAmount();
  amounts.load(node, amount);

  const temperatureMessage = workup.getTemperature();
  if (temperatureMessage) {
    temperature.load(temperatureMessage, node);
  }

  $('.workup_keep_phase', node).text(workup.getKeepPhase());

  const stirringMessage = workup.getStirring();
  if (stirringMessage) {
    stirring.load(stirringMessage, node);
  }

  if (workup.hasTargetPh()) {
    $('.workup_target_ph', node).text(workup.getTargetPh());
  }
  utils.setOptionalBool(
      $('.workup_automated', node),
      workup.hasIsAutomated() ? workup.getIsAutomated() : null);
}

/**
 * Loads a measurement into the given node in a workup.
 * @param {!jQuery} workupNode The div corresponding to the workup whose fields
 *     should be updated.
 * @param {!TemperatureConditions.Measurement} measurement
 */
function loadMeasurement(workupNode, measurement) {
  const node = addMeasurement(workupNode);
  utils.setSelector(
      $('.workup_temperature_measurement_type', node), measurement.getType());
  $('.workup_temperature_measurement_details', node)
      .text(measurement.getDetails());
  const time = measurement.getTime();
  if (time) {
    utils.writeMetric('.workup_temperature_measurement_time', time, node);
  }
  const temperature = measurement.getTemperature();
  if (temperature) {
    utils.writeMetric(
        '.workup_temperature_measurement_temperature', temperature, node);
  }
}

/**
 * Fetches a list of workups defined in the form.
 * @return {!Array<!ReactionWorkup>} workups
 */
function unload() {
  const workups = [];
  $('.workup').each(function(index, node) {
    node = $(node);
    if (!utils.isTemplateOrUndoBuffer(node)) {
      const workup = unloadWorkup(node);
      if (!utils.isEmptyMessage(workup)) {
        workups.push(workup);
      }
    }
  });
  return workups;
}

/**
 * Fetches a single workup from the form.
 * @param {!jQuery} node The div corresponding to the workup to fetch.
 * @return {!ReactionWorkup}
 */
function unloadWorkup(node) {
  const workup = new ReactionWorkup();
  const workupType = utils.getSelectorText($('.workup_type', node)[0]);
  workup.setType(WorkupType[workupType]);
  workup.setDetails(asserts.assertString($('.workup_details', node).text()));

  const duration = utils.readMetric('.workup_duration', new Time(), node);
  if (!utils.isEmptyMessage(duration)) {
    workup.setDuration(duration);
  }

  const input = inputs.unloadInputUnnamed(node);
  if (!utils.isEmptyMessage(input)) {
    workup.setInput(input);
  }

  const amount = amounts.unload($('.workup_amount', node));
  if (!utils.isEmptyMessage(amount)) {
    workup.setAmount(amount);
  }

  const temperatureMessage = temperature.unload(node);
  if (!utils.isEmptyMessage(temperatureMessage)) {
    workup.setTemperature(temperatureMessage);
  }

  workup.setKeepPhase(
      asserts.assertString($('.workup_keep_phase', node).text()));

  const stirringMessage = stirring.unload(node);
  if (!utils.isEmptyMessage(stirringMessage)) {
    workup.setStirring(stirringMessage);
  }

  const targetPh = parseFloat($('.workup_target_ph', node).text());
  if (!isNaN(targetPh)) {
    workup.setTargetPh(targetPh);
  }
  const isAutomated = utils.getOptionalBool($('.workup_automated', node));
  if (isAutomated !== null) {
    workup.setIsAutomated(isAutomated);
  }
  return workup;
}

/**
 * Fetches a single workup temperature measurement from the form.
 * @param {!jQuery} node The div corresponding to the measurement to fetch.
 * @return {!TemperatureConditions.Measurement}
 */
function unloadMeasurement(node) {
  const measurement = new Measurement();
  const measurementType =
      utils.getSelectorText($('.workup_temperature_measurement_type', node)[0]);
  measurement.setType(MeasurementType[measurementType]);
  measurement.setDetails(asserts.assertString(
      $('.workup_temperature_measurement_details', node).text()));
  const time = utils.readMetric(
      '.workup_temperature_measurement_time', new Time(), node);
  if (!utils.isEmptyMessage(time)) {
    measurement.setTime(time);
  }
  const temperature = utils.readMetric(
      '.workup_temperature_measurement_temperature', new Temperature(), node);
  if (!utils.isEmptyMessage(temperature)) {
    measurement.setTemperature(temperature);
  }
  return measurement;
}

/**
 * Adds a new reaction workup section to the form.
 * @return {!jQuery} The newly added parent node for the reaction workup.
 */
function add() {
  const workupNode = utils.addSlowly('#workup_template', $('#workups'));
  const inputNode = $('.workup_input', workupNode);
  // The template for ReactionWorkup.input is taken from Reaction.inputs.
  const workupInputNode = inputs.add(inputNode, ['workup_input']);
  $('.input_addition_order_row', workupInputNode).hide();
  $('.input_addition_time_row', workupInputNode).hide();
  $('.input_addition_duration_row', workupInputNode).hide();
  $('.input_addition_temperature_row', workupInputNode).hide();
  // Adjust heading sizes. Start with the smallest so we don't adjust more than
  // once.
  // TODO(kearnes): This does not affect input components added later.
  $('.h5', workupInputNode).addClass('h6').removeClass('h5');
  $('.h4', workupInputNode).addClass('h5').removeClass('h4');
  $('.h3', workupInputNode).addClass('h4').removeClass('h3');
  // Unlike Reaction.inputs, this ReactionInput has no name.
  $('.input_name', inputNode).hide();
  // Unlike Reaction.inputs, this ReactionInput is not repeated.
  $('.remove', inputNode).hide();

  amounts.init(workupNode);

  // Add live validation handling.
  utils.addChangeHandler(workupNode, () => { validateWorkup(workupNode); });

  // Show/hide fields based on the measurement type.
  const workupTypeSelector = $('.workup_type select', workupNode);
  workupTypeSelector.on('change', function() {
    const workupType = this.options[this.selectedIndex].text;
    if (workupType === 'UNSPECIFIED') {
      $('.workup_keep_phase_row', workupNode).hide();
      $('.workup_target_ph_row', workupNode).hide();
      $('.workup_duration_row', workupNode).hide();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).hide();
      $('.temperature_conditions', workupNode).hide();
      $('.stirring_conditions', workupNode).hide();
    } else if (workupType === 'CUSTOM') {
      $('.workup_keep_phase_row', workupNode).show();
      $('.workup_target_ph_row', workupNode).show();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).show();
      $('.workup_input', workupNode).show();
      $('.temperature_conditions', workupNode).show();
      $('.stirring_conditions', workupNode).show();
    } else if (workupType === 'ADDITION') {
      $('.workup_keep_phase_row', workupNode).hide();
      $('.workup_target_ph_row', workupNode).show();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).show();
      $('.temperature_conditions', workupNode).hide();
      $('.stirring_conditions', workupNode).show();
    } else if (workupType === 'ALIQUOT') {
      $('.workup_keep_phase_row', workupNode).hide();
      $('.workup_target_ph_row', workupNode).hide();
      $('.workup_duration_row', workupNode).hide();
      $('.workup_amount', workupNode).show();
      $('.workup_input', workupNode).hide();
      $('.temperature_conditions', workupNode).hide();
      $('.stirring_conditions', workupNode).hide();
    } else if (workupType === 'TEMPERATURE') {
      $('.workup_keep_phase_row', workupNode).hide();
      $('.workup_target_ph_row', workupNode).hide();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).hide();
      $('.temperature_conditions', workupNode).show();
      $('.stirring_conditions', workupNode).show();
    } else if (workupType === 'CONCENTRATION') {
      $('.workup_keep_phase_row', workupNode).hide();
      $('.workup_target_ph_row', workupNode).hide();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).hide();
      $('.temperature_conditions', workupNode).hide();
      $('.stirring_conditions', workupNode).show();
    } else if (workupType === 'EXTRACTION') {
      $('.workup_keep_phase_row', workupNode).show();
      $('.workup_target_ph_row', workupNode).hide();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).hide();
      $('.temperature_conditions', workupNode).hide();
      $('.stirring_conditions', workupNode).show();
    } else if (workupType === 'FILTRATION') {
      $('.workup_keep_phase_row', workupNode).show();
      $('.workup_target_ph_row', workupNode).hide();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).hide();
      $('.temperature_conditions', workupNode).hide();
      $('.stirring_conditions', workupNode).show();
    } else if (workupType === 'WASH') {
      $('.workup_keep_phase_row', workupNode).hide();
      $('.workup_target_ph_row', workupNode).hide();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).show();
      $('.temperature_conditions', workupNode).hide();
      $('.stirring_conditions', workupNode).show();
    } else if (workupType === 'DRY_IN_VACUUM') {
      $('.workup_keep_phase_row', workupNode).hide();
      $('.workup_target_ph_row', workupNode).hide();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).hide();
      $('.temperature_conditions', workupNode).show();
      $('.stirring_conditions', workupNode).show();
    } else if (workupType === 'DRY_WITH_MATERIAL') {
      $('.workup_keep_phase_row', workupNode).hide();
      $('.workup_target_ph_row', workupNode).hide();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).show();
      $('.temperature_conditions', workupNode).show();
      $('.stirring_conditions', workupNode).show();
    } else if (workupType === 'FLASH_CHROMATOGRAPHY') {
      $('.workup_keep_phase_row', workupNode).hide();
      $('.workup_target_ph_row', workupNode).hide();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).hide();
      $('.temperature_conditions', workupNode).hide();
      $('.stirring_conditions', workupNode).hide();
    } else if (workupType === 'OTHER_CHROMATOGRAPHY') {
      $('.workup_keep_phase_row', workupNode).hide();
      $('.workup_target_ph_row', workupNode).hide();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).hide();
      $('.temperature_conditions', workupNode).hide();
      $('.stirring_conditions', workupNode).hide();
    } else if (workupType === 'SCAVENGING') {
      $('.workup_keep_phase_row', workupNode).hide();
      $('.workup_target_ph_row', workupNode).hide();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).show();
      $('.temperature_conditions', workupNode).hide();
      $('.stirring_conditions', workupNode).show();
    } else if (workupType === 'WAIT') {
      $('.workup_keep_phase_row', workupNode).hide();
      $('.workup_target_ph_row', workupNode).hide();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).hide();
      $('.temperature_conditions', workupNode).hide();
      $('.stirring_conditions', workupNode).hide();
    } else if (workupType === 'STIRRING') {
      $('.workup_keep_phase_row', workupNode).hide();
      $('.workup_target_ph_row', workupNode).hide();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).hide();
      $('.temperature_conditions', workupNode).hide();
      $('.stirring_conditions', workupNode).show();
    } else if (workupType === 'PH_ADJUST') {
      $('.workup_keep_phase_row', workupNode).hide();
      $('.workup_target_ph_row', workupNode).show();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).show();
      $('.temperature_conditions', workupNode).hide();
      $('.stirring_conditions', workupNode).show();
    } else if (workupType === 'DISSOLUTION') {
      $('.workup_keep_phase_row', workupNode).hide();
      $('.workup_target_ph_row', workupNode).hide();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).show();
      $('.temperature_conditions', workupNode).hide();
      $('.stirring_conditions', workupNode).show();
    } else if (workupType === 'DISTILLATION') {
      $('.workup_keep_phase_row', workupNode).hide();
      $('.workup_target_ph_row', workupNode).hide();
      $('.workup_duration_row', workupNode).show();
      $('.workup_amount', workupNode).hide();
      $('.workup_input', workupNode).hide();
      $('.temperature_conditions', workupNode).show();
      $('.stirring_conditions', workupNode).show();
    }
  });
  workupTypeSelector.trigger('change');

  return workupNode;
}

/**
 * Adds a new measurement section to the current workup in the form.
 * @param {!jQuery} node The workup div where the new measurement should be
 *     added.
 * @return {!jQuery} The node of the new measurement div.
 */
function addMeasurement(node) {
  return utils.addSlowly(
      '#workup_temperature_measurement_template',
      $('.workup_temperature_measurements', node));
}

/**
 * Validates a workup as defined in the form.
 * @param {!jQuery} node The div containing to the workup in the form.
 * @param {?jQuery=} validateNode The target div for validation results.
 */
function validateWorkup(node, validateNode = null) {
  const workup = unloadWorkup(node);
  utils.validate(workup, 'ReactionWorkup', node, validateNode);
}
