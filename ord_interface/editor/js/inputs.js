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

goog.module('ord.inputs');
goog.module.declareLegacyNamespace();

const asserts = goog.require('goog.asserts');

const JspbMap = goog.requireType('jspb.Map');

const compounds = goog.require('ord.compounds');
const crudes = goog.require('ord.crudes');
const utils = goog.require('ord.utils');

const FlowRate = goog.require('proto.ord.FlowRate');
const ReactionInput = goog.require('proto.ord.ReactionInput');
const AdditionDevice = goog.require('proto.ord.ReactionInput.AdditionDevice');
const AdditionDeviceType = goog.require('proto.ord.ReactionInput.AdditionDevice.AdditionDeviceType');
const AdditionSpeed = goog.require('proto.ord.ReactionInput.AdditionSpeed');
const AdditionSpeedType = goog.require('proto.ord.ReactionInput.AdditionSpeed.AdditionSpeedType');
const Temperature = goog.require('proto.ord.Temperature');
const Time = goog.require('proto.ord.Time');

exports = {
  load,
  loadInputUnnamed,
  unload,
  unloadInputUnnamed,
  add,
  validateInput,
  addInputByString,
};


/**
 * Adds and populates the reaction input sections in the form.
 * @param {!JspbMap<string, !ReactionInput>} inputs
 */
function load(inputs) {
  inputs.forEach(function(input, name) {
    loadInput($('#inputs'), name, input);
  });
  utils.updateSidebar();
}

/**
 * Adds and populates a single reaction input section in the form.
 * @param {!jQuery} root Root node for the reaction input.
 * @param {string} name The name of this input.
 * @param {!ReactionInput} input
 * @return {!jQuery} The new input node.
 */
function loadInput(root, name, input) {
  const node = add(root);
  loadInputUnnamed(node, input);
  $('.input_name', node).text(name);
  return node;
}

/**
 * Adds and populates a single reaction input section in the form according to
 * a short-hand string describing a stock solution.
 * @param {!jQuery} root Root node for the reaction input.
 */
function addInputByString(root) {
  const string =
      prompt(`Please describe the input in one of the following forms:
(1) [AMOUNT] of [NAME]
(2) [AMOUNT] of [CONCENTRATION] [SOLUTE] in [SOLVENT]`);
  if (!(string)) {
    return;
  }
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/resolve/input');
  xhr.responseType = 'arraybuffer';
  xhr.onload = () => {
    if (xhr.status === 409) {
      const decoder = new TextDecoder('utf-8');
      asserts.assertInstanceof(xhr.response, ArrayBuffer);  // Type hint.
      alert('Could not parse: ' + decoder.decode(xhr.response));
    } else {
      asserts.assertInstanceof(xhr.response, ArrayBuffer);  // Type hint.
      const bytes = new Uint8Array(xhr.response);
      const input = ReactionInput.deserializeBinary(bytes);
      if (input) {
        const input_node = loadInput(root, asserts.assertString(string), input);
        $('.edittext', input_node).trigger('blur');
      }
    }
  };
  xhr.send(string);
}

/**
 * Adds and populates a single reaction input section in the form without
 * assigning it a name.
 * @param {!jQuery} node Root node for the reaction input.
 * @param {?ReactionInput} input
 * @return {!jQuery} The original root node.
 */
function loadInputUnnamed(node, input) {
  if (!input) {
    return node;
  }
  const componentsList = input.getComponentsList();
  compounds.load(node, componentsList);

  const crudeComponentsList = input.getCrudeComponentsList();
  crudes.load(node, crudeComponentsList);

  const additionOrder = input.getAdditionOrder();
  if (additionOrder !== 0) {
    $('.input_addition_order', node).text(additionOrder);
  }

  const additionTime = input.getAdditionTime();
  if (additionTime) {
    utils.writeMetric('.input_addition_time', additionTime, node);
  }
  const additionSpeed = input.getAdditionSpeed();
  if (additionSpeed) {
    utils.setSelector(
        $('.input_addition_speed_type', node), additionSpeed.getType());
    $('.input_addition_speed_details', node).text(additionSpeed.getDetails());
  }
  const additionDevice = input.getAdditionDevice();
  if (additionDevice) {
    utils.setSelector(
        $('.input_addition_device_type', node), additionDevice.getType());
    $('.input_addition_device_details', node).text(additionDevice.getDetails());
  }
  const duration = input.getAdditionDuration();
  if (duration) {
    utils.writeMetric('.input_addition_duration', duration, node);
  }
  const temperature = input.getAdditionTemperature();
  if (temperature) {
    utils.writeMetric('.input_addition_temperature', temperature, node);
  }
  const flowRate = input.getFlowRate();
  if (flowRate) {
    utils.writeMetric('.input_flow_rate', flowRate, node);
  }
  return node;
}

/**
 * Fetches the reaction inputs defined in the form and adds them to `inputs`.
 * @param {!JspbMap<string, !ReactionInput>} inputs
 */
function unload(inputs) {
  $('#inputs > div.input').each(function(index, node) {
    node = $(node);
    if (!utils.isTemplateOrUndoBuffer(node)) {
      unloadInput(inputs, node);
    }
  });
}

/**
 * Fetches a single reaction input defined in the form and adds it to `inputs`.
 * @param {!JspbMap<string, !ReactionInput>} inputs
 * @param {!jQuery} node Root node for the reaction input.
 */
function unloadInput(inputs, node) {
  const name = $('.input_name', node).text();
  const input = unloadInputUnnamed(node);
  if (name || !utils.isEmptyMessage(input)) {
    inputs.set(asserts.assertString(name), input);
  }
}

/**
 * Fetches a single reaction input defined in the form.
 * @param {!jQuery} node Root node for the reaction input.
 * @return {!ReactionInput}
 */
function unloadInputUnnamed(node) {
  const input = new ReactionInput();

  const componentsList = compounds.unload(node);
  if (componentsList.some(e => !utils.isEmptyMessage(e))) {
    input.setComponentsList(componentsList);
  }

  const crudeComponentsList = crudes.unload(node);
  if (crudeComponentsList.some(e => !utils.isEmptyMessage(e))) {
    input.setCrudeComponentsList(crudeComponentsList);
  }

  const additionOrder = parseInt($('.input_addition_order', node).text(), 10);
  if (!isNaN(additionOrder)) {
    input.setAdditionOrder(additionOrder);
  }
  const additionTime =
      utils.readMetric('.input_addition_time', new Time(), node);
  if (!utils.isEmptyMessage(additionTime)) {
    input.setAdditionTime(additionTime);
  }

  const additionSpeed = new AdditionSpeed();
  const additionSpeedType =
      utils.getSelectorText($('.input_addition_speed_type', node)[0]);
  additionSpeed.setType(AdditionSpeedType[additionSpeedType]);
  additionSpeed.setDetails(
      asserts.assertString($('.input_addition_speed_details', node).text()));
  if (!utils.isEmptyMessage(additionSpeed)) {
    input.setAdditionSpeed(additionSpeed);
  }

  const additionDevice = new AdditionDevice();
  const additionDeviceType =
      utils.getSelectorText($('.input_addition_device_type', node)[0]);
  additionDevice.setType(AdditionDeviceType[additionDeviceType]);
  additionDevice.setDetails(
      asserts.assertString($('.input_addition_device_details', node).text()));
  if (!utils.isEmptyMessage(additionDevice)) {
    input.setAdditionDevice(additionDevice);
  }

  const additionDuration =
      utils.readMetric('.input_addition_duration', new Time(), node);
  if (!utils.isEmptyMessage(additionDuration)) {
    input.setAdditionDuration(additionDuration);
  }

  const additionTemperature =
      utils.readMetric('.input_addition_temperature', new Temperature(), node);
  if (!utils.isEmptyMessage(additionTemperature)) {
    input.setAdditionTemperature(additionTemperature);
  }

  const flowRate = utils.readMetric('.input_flow_rate', new FlowRate(), node);
  if (!utils.isEmptyMessage(flowRate)) {
    input.setFlowRate(flowRate);
  }

  return input;
}

/**
 * Adds a reaction input section to the form.
 * @param {!jQuery} root Parent node for reaction inputs.
 * @param {?Array<string>=} classes Additional classes to add to the new node.
 * @return {!jQuery} The newly added parent node for the reaction input.
 */
function add(root, classes = null) {
  const node = utils.addSlowly('#input_template', $(root));
  if (Array.isArray(classes) && classes.length) {
    node.addClass(classes);
  }
  utils.updateSidebar();
  // Add live validation handling.
  utils.addChangeHandler(node, () => {
    validateInput(node);
  });
  // Update the sidebar when the input name is changed.
  const nameNode = node.find('.input_name').first();
  nameNode.on('blur', utils.updateSidebar);
  return node;
}

/**
 * Validates a reaction input defined in the form.
 * @param {!jQuery} node Root node for the reaction input.
 * @param {?jQuery=} validateNode Target node for validation results.
 */
function validateInput(node, validateNode = null) {
  const input = unloadInputUnnamed(node);
  utils.validate(input, 'ReactionInput', node, validateNode);
}
