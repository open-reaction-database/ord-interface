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

goog.module('ord.utils');
goog.module.declareLegacyNamespace();

const Message = goog.require('jspb.Message');

const asserts = goog.require('goog.asserts');

/** @suppress {extraRequire} */
const enums = goog.require('ord.enums');  // Used by nameToProto.

const Dataset = goog.require('proto.ord.Dataset');
const Current = goog.require('proto.ord.Current');
const FlowRate = goog.require('proto.ord.FlowRate');
const Length = goog.require('proto.ord.Length');
const Mass = goog.require('proto.ord.Mass');
const Moles = goog.require('proto.ord.Moles');
const Percentage = goog.require('proto.ord.Percentage');
const Pressure = goog.require('proto.ord.Pressure');
const Reaction = goog.require('proto.ord.Reaction');
const Temperature = goog.require('proto.ord.Temperature');
const Time = goog.require('proto.ord.Time');
const Voltage = goog.require('proto.ord.Voltage');
const Volume = goog.require('proto.ord.Volume');
const Wavelength = goog.require('proto.ord.Wavelength');

/**
 * Messages with 'units' fields.
 * @typedef {!Current|
 *           !FlowRate|
 *           !Length|
 *           !Mass|
 *           !Moles|
 *           !Pressure|
 *           !Temperature|
 *           !Time|
 *           !Voltage|
 *           !Volume|
 *           !Wavelength}
 */
let UnitMessage;

exports = {
  addChangeHandler,
  addSlowly,
  clean,
  compareDataset,
  getDataset,
  getOptionalBool,
  getReactionById,
  getSelectorText,
  freeze,
  initOptionalBool,
  initSelector,
  initValidateNode,
  isEmptyMessage,
  isTemplateOrUndoBuffer,
  listen,
  prepareFloat,
  putDataset,
  ready,
  readMetric,
  removeSlowly,
  setupObserver,
  setOptionalBool,
  setSelector,
  setTextFromFile,
  showOptionalSection,
  toggleAutosave,
  toggleValidateMessage,
  undoSlowly,
  updateSidebar,
  validate,
  writeMetric,
};

// Remember the dataset and reaction we are editing.
const session = {
  fileName: null,
  dataset: null,
  index: null,             // Ordinal position of the Reaction in its Dataset.
  observer: null,          // IntersectionObserver used for the sidebar.
  navSelectors: {},        // Dictionary from navigation to section.
  timers: {'short': null}  // A timer used by autosave.
};
// Export session, because it's used by test.js.
exports.session = session;

const FLOAT_PATTERN = /^-?(?:\d+|\d+\.\d*|\d*\.\d+)(?:[eE]-?\d+)?$/;
const INTEGER_PATTERN = /^-?\d+$/;

/**
 * Sets the `ready` value to true.
 * @suppress {undefinedVars}
 */
function ready() {
  $('body').attr('ready', 'true');
}

/**
 * Shows the 'save' button.
 */
function dirty() {
  $('#save').css('visibility', 'visible');
}

/**
 * Hides the 'save' button.
 */
function clean() {
  const matcher = $('#save');
  matcher.css('visibility', 'hidden');
  matcher.text('save');
}

/**
 * Adds a change handler to the given node that shows the 'save' button when
 * the node text is edited.
 * @param {!jQuery} node
 */
function listen(node) {
  addChangeHandler(node, dirty);
  $('.edittext', node).on('focus', event => selectText(event.target));
  $('.floattext', node).on('blur', event => checkFloat(event.target));
  $('.integertext', node).on('blur', event => checkInteger(event.target));
}

/**
 * Clicks the 'save' button if ready for a save.
 */
function clickSave() {
  // Only save if there are unsaved changes still to be saved -- hence save
  // button visible -- and if ready for a save (not in the process of saving
  // already).
  const saveButton = $('#save');
  if (saveButton.css('visibility') === 'visible' &&
      saveButton.text() === 'save') {
    saveButton.trigger('click');
    // Validate on autosave.
    $('#reaction_validate_button').trigger('click');
  }
}

/**
 * Toggles autosave being active.
 */
function toggleAutosave() {
  const matcher = $('#toggle_autosave');
  // We keep track of timers by holding references, only if they're active.
  if (!session.timers.short) {
    // Enable a simple timer that saves periodically.
    session.timers.short =
        setInterval(clickSave, 1000 * 15);  // Save after 15 seconds
    matcher.text('autosave: on');
  } else {
    // Stop the interval timer itself, then remove reference in order to
    // properly later detect that it's stopped.
    clearInterval(session.timers.short);
    session.timers.short = null;
    matcher.text('autosave: off');
  }
}

/**
 * Adds an instance of `template` to the root node.
 * @param {string} template A jQuery selector.
 * @param {!jQuery} root A jQuery object.
 * @return {!jQuery} The new copy of the template.
 *
 * NOTE(kearnes): .tooltip() is part of jQuery UI.
 * @suppress {missingProperties}
 */
function addSlowly(template, root) {
  const node = $(template).clone();
  node.removeAttr('id');
  root.append(node);
  node.show('slow');
  dirty();
  listen(node);
  $('[data-toggle="tooltip"]', node).tooltip();
  return node;
}

/**
 * Removes from the DOM the nearest ancestor element matching the pattern.
 * @param {!jQuery} button The element from which to start the search.
 * @param {string} pattern The pattern for the element to remove.
 */
function removeSlowly(button, pattern) {
  const node = $(button.closest(pattern));
  // Must call necessary validators only after the node is removed,
  // but we can only figure out which validators these are before removal.
  // We do so, and after removal, click the corresponding buttons to trigger
  // validation.
  let buttonsToClick = $();
  node.parents('fieldset').each(function() {
    buttonsToClick =
        buttonsToClick.add($(this).children('legend').find('.validate_button'));
  });
  makeUndoable(node);
  node.hide('slow', function() {
    buttonsToClick.trigger('click');
    updateSidebar();
  });
  dirty();
}

/**
 * Reverses the hide() in the most recent invocation of removeSlowly().
 * Removes the node's "undo" button. Does not trigger validation.
 */
function undoSlowly() {
  $('.undoable').removeClass('undoable').show('slow');
  $('.undo').not('#undo_template').hide('slow', function() {
    $(this).remove();
    updateSidebar();
  });
  dirty();
}

/**
 * Marks the given node for possible future undo. Adds an "undo" button to do
 * it. Deletes any preexisting undoable nodes and undo buttons.
 * @param {!jQuery} node The DOM fragment to hide and re-show.
 */
function makeUndoable(node) {
  $('.undoable').remove();
  node.addClass('undoable');
  $('.undo').not('#undo_template').remove();
  const button = $('#undo_template').clone();
  button.removeAttr('id');
  node.after(button);
  button.show('slow');
}

/**
 * On button click, hides the button and displays an optional input section.
 * @param {!jQuery} button The button that was clicked.
 * @param {string} target ID of the target element.
 */
function showOptionalSection(button, target) {
  button.hide();
  $('#' + target).toggle('slow');
}

/**
 * Adds and populates a <select/> node according to its data-proto type
 * declaration.
 * @param {!jQuery} node A node containing a `data-proto` attribute.
 */
function initSelector(node) {
  const protoName = node.attr('data-proto');
  const protoEnum = nameToProto(asserts.assertString(protoName));
  if (!protoEnum) {
    console.log('missing require: "' + protoName + '"');
  }
  const types = Object.entries(asserts.assertObject(protoEnum));
  const select = $('<select>');
  for (let i = 0; i < types.length; i++) {
    const option = $('<option>').text(types[i][0]);
    option.attr('value', types[i][1]);
    if (types[i][0] === 'UNSPECIFIED') {
      option.attr('selected', 'selected');
    }
    select.append(option);
  }
  node.append(select);
}

/**
 * Sets up a three-way popup (true/false/unspecified).
 * @param {!jQuery} node Target node for the new <select/> element.
 */
function initOptionalBool(node) {
  const select = $('<select>');
  const options = ['UNSPECIFIED', 'TRUE', 'FALSE'];
  for (let i = 0; i < options.length; i++) {
    const option = $('<option>').text(options[i]);
    option.attr('value', options[i]);
    if (options[i] === 'UNSPECIFIED') {
      option.attr('selected', 'selected');
    }
    select.append(option);
  }
  node.append(select);
}

/**
 * Sets up a validator div (button, status indicator, error list, etc.) by
 * inserting contents into a div in reaction.html.
 * @param {!jQuery} oldNode Target node for the new validation elements.
 */
function initValidateNode(oldNode) {
  const newNode = $('#validate_template').clone();
  oldNode.append(newNode.children());
}

/**
 * Fetches a reaction as a serialized Reaction proto.
 * @param {string} reactionId The ID of the Reaction to fetch.
 * @return {!Promise<!Uint8Array>}
 */
function getReactionById(reactionId) {
  return new Promise(resolve => {
    const xhr = new XMLHttpRequest();
    xhr.open('GET', '/reaction/id/' + reactionId + '/proto');
    xhr.responseType = 'arraybuffer';
    xhr.onload = function() {
      asserts.assertInstanceof(xhr.response, ArrayBuffer);  // Type hint.
      const bytes = new Uint8Array(xhr.response);
      const reaction = Reaction.deserializeBinary(bytes);
      resolve(reaction);
    };
    xhr.send();
  });
}

/**
 * Converts a Message_Field name from a data-proto attribute into a proto class.
 * @param {string} protoName Underscore-delimited protocol buffer field name,
 *     such as Reaction_provenance.
 * @return {?typeof Message}
 */
function nameToProto(protoName) {
  let clazz = proto.ord;
  protoName.split('_').forEach(function(name) {
    clazz = clazz[name];
    if (!clazz) {
      return null;
    }
  });
  return clazz;
}

/**
 * Selects the contents of the given node.
 * @param {!Node} node
 */
function selectText(node) {
  const range = document.createRange();
  range.selectNodeContents(node);
  const selection = window.getSelection();
  selection.removeAllRanges();
  selection.addRange(range);
}

/**
 * Determines if the text entered in a float input is valid by detecting any
 * characters besides 0-9, a single period to signify a decimal, and a
 * leading hyphen. Also supports scientific notation with either 'e' or 'E'.
 * @param {!jQuery} node
 */
function checkFloat(node) {
  const stringValue = $(node).text().trim();
  if (stringValue === '') {
    $(node).removeClass('invalid');
  } else if (FLOAT_PATTERN.test(stringValue)) {
    $(node).removeClass('invalid');
  } else {
    $(node).addClass('invalid');
  }
}

/**
 * Determines if the text entered in an integer input is valid by forbidding
 * any characters besides 0-9 and a leading hyphen.
 * @param {!jQuery} node
 */
function checkInteger(node) {
  const stringValue = $(node).text().trim();
  if (stringValue === '') {
    $(node).removeClass('invalid');
  } else if (INTEGER_PATTERN.test(stringValue)) {
    $(node).removeClass('invalid');
  } else {
    $(node).addClass('invalid');
  }
}

/**
 * Prepares a floating point value for display in the form.
 * @param {number} value
 * @return {number}
 */
function prepareFloat(value) {
  // Round to N significant digits; this avoid floating point precision issues
  // that can be quite jarring to users.
  //
  // See:
  //   * https://stackoverflow.com/a/3644302
  //   * https://medium.com/swlh/ed74c471c1b8
  //   * https://stackoverflow.com/a/19623253
  const precision = 7;
  return parseFloat(value.toPrecision(precision));
}

/**
 * Adds a jQuery handler to a node such that the handler is run once whenever
 * data entry within that node is changed, *except through remove* -- this must
 * be handled manually. (This prevents inconsistent timing in the ordering of
 * the element being removed and the handler being called.)
 * @param {!jQuery} node
 * @param {!Function} handler
 */
function addChangeHandler(node, handler) {
  // For textboxes
  node.on('blur', '.edittext', handler);
  // For selectors, optional bool selectors,
  // and checkboxes/radio buttons/file upload, respectively
  node.on('change', '.selector, .optional_bool, input', handler);
  // For add buttons
  node.on('click', '.add', handler);
  // For validate divs.
  node.on('click', '.validate', handler);
}

/**
 * Generic validator for many message types, not just reaction.
 * NOTE: This function does not commit or save anything!
 * @param {!Message} message The proto to validate.
 * @param {string} messageTypeString The message type.
 * @param {!jQuery} node Parent node for the unloaded message.
 * @param {?jQuery} validateNode Target div for validation output.
 *
 * NOTE(kearnes): serializeBinary is not defined in the base class.
 * TODO(kearnes): Annotate `errors` and `warnings` properties on response.
 * @suppress {missingProperties}
 */
function validate(message, messageTypeString, node, validateNode) {
  // eg message is a type of reaction, messageTypeString = 'Reaction'
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/dataset/proto/validate/' + messageTypeString);
  const binary = message.serializeBinary();
  if (!validateNode) {
    validateNode = $('.validate', node).first();
  }
  xhr.responseType = 'json';
  xhr.onload = function() {
    const validationOutput = xhr.response;
    const errors = validationOutput.errors;
    const warnings = validationOutput.warnings;
    // Add client-side validation errors.
    node.find('.invalid').each(function() {
      const invalidName = $(this).attr('class').split(' ')[0];
      errors.push('Value for ' + invalidName + ' is invalid');
    });
    const statusNode = $('.validate_status', validateNode);
    const messageNode = $('.validate_message', validateNode);
    if (errors.length) {
      statusNode.show();
      statusNode.text(' ' + errors.length);
      messageNode.show();
      messageNode.html('<ul></ul>');
      for (let index = 0; index < errors.length; index++) {
        const error = errors[index];
        const errorNode = $('<li></li>');
        errorNode.text(error);
        $('ul', messageNode).append(errorNode);
      }
    } else {
      statusNode.hide();
      messageNode.html('');
      messageNode.hide();
    }
    const warningStatusNode = $('.validate_warning_status', validateNode);
    const warningMessageNode = $('.validate_warning_message', validateNode);
    if (warnings.length) {
      warningStatusNode.show();
      warningStatusNode.text(' ' + warnings.length);
      warningMessageNode.show();
      warningMessageNode.html('<ul></ul>');
      for (let index = 0; index < warnings.length; index++) {
        const warning = warnings[index];
        const warningNode = $('<li></li>');
        warningNode.text(warning);
        $('ul', warningMessageNode).append(warningNode);
      }
    } else {
      warningStatusNode.hide();
      warningMessageNode.html('');
      warningMessageNode.hide();
    }
  };
  xhr.send(binary);
}

/**
 * Toggles the visibility of the 'validate' button for a given node.
 * @param {!jQuery} target
 */
function toggleValidateMessage(target) {
  switch (target.css('visibility')) {
    case 'visible':
      target.css('visibility', 'hidden');
      break;
    case 'hidden':
      target.css('visibility', 'visible');
      break;
  }
}

/**
 * Downloads a dataset as a serialized Dataset proto.
 * @param {string} fileName The name of the dataset to fetch.
 * @return {!Promise<!Uint8Array>}
 */
function getDataset(fileName) {
  return new Promise(resolve => {
    const xhr = new XMLHttpRequest();
    xhr.open('GET', '/dataset/proto/read/' + fileName);
    xhr.responseType = 'arraybuffer';
    xhr.onload = function() {
      asserts.assertInstanceof(xhr.response, ArrayBuffer);  // Type hint.
      const bytes = new Uint8Array(xhr.response);
      const dataset = Dataset.deserializeBinary(bytes);
      resolve(dataset);
    };
    xhr.send();
  });
}

/**
 * Uploads a serialized Dataset proto.
 * @param {string} fileName The name of the new dataset.
 * @param {!Dataset} dataset
 */
function putDataset(fileName, dataset) {
  $('#save').text('saving');
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/dataset/proto/write/' + fileName);
  const binary = dataset.serializeBinary();
  xhr.onload = clean;
  xhr.send(binary);
}

/**
 * Compares a local Dataset to a Dataset on the server (used for testing).
 * @param {string} fileName The name of a dataset on the server.
 * @param {!Dataset} dataset A local Dataset.
 * @return {!Promise<!Uint8Array>}
 */
async function compareDataset(fileName, dataset) {
  return new Promise((resolve, reject) => {
    const xhr = new XMLHttpRequest();
    xhr.open('POST', '/dataset/proto/compare/' + fileName);
    const binary = dataset.serializeBinary();
    xhr.onload = () => {
      if (xhr.status === 200) {
        resolve();
      } else {
        reject();
      }
    };
    xhr.onerror = reject;
    xhr.send(binary);
  });
}

/**
 * Checks if the argument represents an empty protobuf message (that is, the
 * argument's nested arrays only contains null or empty values), or is null or
 * undefined. We use this check on both primitives and arrays/messages.
 * NOTE: Unlike other primitive types, using a setter to set a oneof string
 * field to “” causes the message to include the field and “”, which would be
 * unwanted. So we instead claim that empty strings are empty messages. (Hence
 * we don’t set _any_ empty string.)
 * NOTE: In a submessage, setting a meaningful value (e.g. optional float to 0)
 * will result in a non-null/undefined value in the submessage array. So, if
 * the array of a submessage only contains null and undefined vals, we can
 * assume that the message is truly “empty” (that is, doesn’t have anything
 * meaningful that is set) and can be omitted when constructing the surrounding
 * message.
 * @param {!Message} obj The object to test.
 * @return {boolean} Whether the message is empty.
 *
 * NOTE(kearnes): serializeBinary is not defined in the base class.
 * @suppress {missingProperties}
 */
function isEmptyMessage(obj) {
  const empty = new obj.constructor();
  // Compare binary encodings to cover optional fields.
  return JSON.stringify(obj.serializeBinary()) ===
      JSON.stringify(empty.serializeBinary());
}

/**
 * Supports unload() operations by filtering spurious selector matches due
 * either to DOM templates or elements the user has removed undoably.
 * @param {!jQuery} node The DOM node to test for spuriousness.
 * @return {boolean} True means ignore this node.
 */
function isTemplateOrUndoBuffer(node) {
  return !!(node.attr('id') || node.hasClass('undoable'));
}

/**
 * Unpacks a (value, units, precision) tuple into the given type.
 * @param {string} prefix The prefix for element attributes.
 * @param {TYPE} proto A protocol buffer with `value`, `precision`,
 *     and `units` fields.
 * @param {?jQuery=} node The node containing the tuple.
 * @return {TYPE} The updated protocol buffer. Note that the message
 *     is modified in-place.
 * @template TYPE
 */
function readMetric(prefix, proto, node = null) {
  const value = parseFloat($(prefix + '_value', node).text());
  if (!isNaN(value)) {
    proto.setValue(value);
  }
  if (proto.setUnits) {
    // proto.ord.Percentage doesn't have units.
    const unitsNode = $(prefix + '_units', node);
    const unitsName = unitsNode.attr('data-proto');
    const unitsEnum = nameToProto(asserts.assertString(unitsName));
    const units = getSelectorText(unitsNode[0]);
    proto.setUnits(unitsEnum[units]);
  }
  const precision = parseFloat($(prefix + '_precision', node).text());
  if (!isNaN(precision)) {
    proto.setPrecision(precision);
  }
  return proto;
}

/**
 * Unpacks a (value, units, precision) tuple into form elements.
 * @param {string} prefix The prefix for element attributes.
 * @param {?UnitMessage|?Percentage} proto A protocol buffer with
 *    `value`, `precision`, and `units` fields. (Percentage is allowed
 *    as a special case even though it doesn't have a `units` field.)
 * @param {?jQuery=} node The target node for the tuple.
 */
function writeMetric(prefix, proto, node = null) {
  if (!proto) {
    return;
  }
  if (proto.hasValue()) {
    $(prefix + '_value', node).text(proto.getValue());
  }
  if (proto.getUnits) {
    // proto.ord.Percentage doesn't have units.
    setSelector($(prefix + '_units', node), proto.getUnits());
  }
  if (proto.hasPrecision()) {
    $(prefix + '_precision', node).text(proto.getPrecision());
  }
}

/**
 * Prompts the user to upload a file and sets the target node text with its
 * contents.
 * @param {!jQuery} identifierNode The node to update with the file contents.
 * @param {string} valueClass The class containing `identifierNode`.
 */
function setTextFromFile(identifierNode, valueClass) {
  const input = document.createElement('input');
  input.setAttribute('type', 'file');
  input.onchange = (event => {
    const file = event.target.files[0];
    const reader = new FileReader();
    reader.readAsText(file);
    reader.onload = readerEvent => {
      const contents = readerEvent.target.result;
      $('.' + valueClass, identifierNode).text(contents);
    };
  });
  input.click();
}

/**
 * Selects an <option/> under a <select/>.
 * @param {!jQuery} node A <select/> element.
 * @param {number} value
 */
function setSelector(node, value) {
  $('option', node).first().removeAttr('selected');
  $('option[value=' + value + ']', node).first().attr('selected', 'selected');
}

/**
 * Finds the selected <option/> and returns its text.
 * @param {!Element} node A node containing one or more <select/> elements.
 * @return {string}
 */
function getSelectorText(node) {
  const selectorElement = node.getElementsByTagName('select')[0];
  asserts.assertInstanceof(selectorElement, HTMLSelectElement);  // Type hint.
  return selectorElement.options[selectorElement.selectedIndex].text;
}

/**
 * Sets the value of a three-way popup (true/false/unspecified).
 * @param {!jQuery} node A node containing a three-way selector.
 * @param {boolean|null} value The value to select.
 */
function setOptionalBool(node, value) {
  $('option', node).removeAttr('selected');
  if (value === true) {
    $('option[value=TRUE]', node).attr('selected', 'selected');
  }
  if (value === false) {
    $('option[value=FALSE]', node).attr('selected', 'selected');
  }
  if (value == null) {
    $('option[value=UNSPECIFIED]', node).attr('selected', 'selected');
  }
}

/**
 * Fetches the value of a three-way popup (true/false/unspecified).
 * @param {!jQuery} node A node containing a three-way selector.
 * @return {boolean|null}
 */
function getOptionalBool(node) {
  const value = $('select', node).val();
  if (value === 'TRUE') {
    return true;
  }
  if (value === 'FALSE') {
    return false;
  }
  return null;
}

/**
 * Switches the UI into a read-only mode. This is irreversible.
 */
function freeze() {
  // Hide the header buttons...
  $('#header_buttons').children().hide();
  // ...except for "download".
  $('#download').show();
  $('#identity').hide();
  $('select').attr('disabled', 'true');
  $('input:radio').prop('disabled', 'true');
  $('.validate').hide();
  $('.add').hide();
  $('.remove').hide();
  $('.text_upload').hide();
  $('#provenance_created button').hide();
  $('.edittext').each((i, x) => {
    const node = $(x);
    node.attr('contenteditable', 'false');
    node.css('background-color', '#ebebe4');
  });
}

/**
 * Highlights navigation buttons in the sidebar corresponding to visible
 * sections. Used as a callback function for the IntersectionObserver.
 * @param {!Array<!IntersectionObserverEntry>} entries
 */
function observerCallback(entries) {
  entries.forEach(entry => {
    const target = $(entry.target);
    let section;
    if (target[0].hasAttribute('input_name')) {
      section = target.attr('input_name');
    } else {
      section = target.attr('id').split('_')[1];
    }
    if (entry.isIntersecting) {
      session.navSelectors[section].css('background-color', 'lightblue');
    } else {
      session.navSelectors[section].css('background-color', '');
    }
  });
}

/**
 * Sets up the IntersectionObserver used to highlight navigation buttons
 * in the sidebar.
 */
function setupObserver() {
  const headerSize = $('#header').outerHeight();
  const observerOptions = {rootMargin: '-' + headerSize + 'px 0px 0px 0px'};
  session.observer =
      new IntersectionObserver(observerCallback, observerOptions);
  updateObserver();
}

/**
 * Updates the set of elements watched by the IntersectionObserver.
 */
function updateObserver() {
  if (!session.observer) {
    return;  // Do nothing until setupObserver has been run.
  }
  session.observer.disconnect();
  $('.section:visible').not('.workup_input').each(function() {
    session.observer.observe(this);
  });
  // Index the selector controls.
  session.navSelectors = {};
  $('.navSection').each((index, selector) => {
    selector = $(selector);
    const section = selector.attr('data-section');
    session.navSelectors[section] = selector;
  });
  $('.inputNavSection').each((index, selector) => {
    selector = $(selector);
    const section = selector.attr('input_name');
    session.navSelectors[section] = selector;
  });
}

/**
 * Scrolls the viewport to the selected input.
 * @param {!Event} event
 */
function scrollToInput(event) {
  const section = $(event.target).attr('input_name');
  const target = $('.input[input_name=\'' + section + '\']');
  target[0].scrollIntoView({behavior: 'smooth'});
}

/**
 * Updates the input entries in the sidebar.
 */
function updateSidebar() {
  $('#navInputs').empty();
  $('.input:visible').not('.workup_input').each(function(index) {
    const node = $(this);
    let name = node.find('.input_name').first().text();
    if (name === '') {
      name = '(Input #' + (index + 1) + ')';
    }
    node.attr('input_name', 'INPUT-' + name);
    const navNode = $('<div>&#8226; ' + name + '</div>');
    navNode.addClass('inputNavSection');
    navNode.attr('input_name', 'INPUT-' + name);
    $('#navInputs').append(navNode);
    navNode.on('click', scrollToInput);
  });
  updateObserver();
}
