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

goog.module('ord.data');
goog.module.declareLegacyNamespace();
exports = {
  addData,
  loadData,
  unloadData,
};

const asserts = goog.require('goog.asserts');

const uploads = goog.require('ord.uploads');
const utils = goog.require('ord.utils');

const Data = goog.require('proto.ord.Data');
const KindCase = goog.require('proto.ord.Data.KindCase');

// Freely create radio button groups by generating new input names.
let radioGroupCounter = 0;

/**
 * Adds a new Data section to the form.
 * @param {!jQuery} parentNode Parent node.
 * @return {!jQuery} The newly added node for the Data record.
 */
function addData(parentNode) {
  const target = parentNode.children('fieldset').first();
  const node = utils.addSlowly('#data_template', target);
  const typeButtons = $('input[type=\'radio\']', node);
  typeButtons.attr('name', 'data_' + radioGroupCounter++);
  typeButtons.on('change', function() {
    if ((this.value === 'text') || (this.value === 'number') ||
        (this.value === 'url')) {
      $('.data_text', node).show();
      $('.data_uploader', node).hide();
    } else {
      $('.data_text', node).hide();
      $('.data_uploader', node).show();
    }
    if (this.value === 'number') {
      $('.data_text', node).addClass('floattext');
    } else {
      $('.data_text', node).removeClass('floattext');
    }
  });
  uploads.initialize(node);
  return node;
}

/**
 * Populates an existing Data section in the form.
 * @param {!jQuery} node Root node.
 * @param {?Data} data
 */
function loadData(node, data) {
  if (!data) {
    return;
  }
  $('.data_description', node).text(data.getDescription());
  $('.data_format', node).text(data.getFormat());
  let value;
  switch (data.getKindCase()) {
    case KindCase.FLOAT_VALUE:
      value = data.getFloatValue().toString();
      if (value.indexOf('.') === -1) {
        value = value.concat('.0');
      }
      $('.data_text', node).show();
      $('.data_uploader', node).hide();
      $('.data_text', node).text(value);
      $('input[value=\'number\']', node).prop('checked', true);
      break;
    case KindCase.INTEGER_VALUE:
      value = data.getIntegerValue();
      $('.data_text', node).show();
      $('.data_uploader', node).hide();
      $('.data_text', node).text(value);
      $('input[value=\'number\']', node).prop('checked', true);
      break;
    case KindCase.BYTES_VALUE:
      value = data.getBytesValue();
      $('.data_text', node).hide();
      $('.data_uploader', node).show();
      asserts.assertInstanceof(value, Uint8Array);  // Type hint.
      uploads.load(node, value);
      $('input[value=\'upload\']', node).prop('checked', true);
      break;
    case KindCase.STRING_VALUE:
      value = data.getStringValue();
      $('.data_text', node).show();
      $('.data_uploader', node).hide();
      $('.data_text', node).text(value);
      $('input[value=\'text\']', node).prop('checked', true);
      break;
    case KindCase.URL:
      value = data.getUrl();
      $('.data_text', node).show();
      $('.data_uploader', node).hide();
      $('.data_text', node).text(value);
      $('input[value=\'url\']', node).prop('checked', true);
      break;
    default:
      break;
  }
}

/**
 * Fetches a Data section from the form.
 * @param {!jQuery} node Root node of the Data section to fetch.
 * @return {!Data}
 */
function unloadData(node) {
  const data = new Data();
  const description = $('.data_description', node).text();
  data.setDescription(asserts.assertString(description));
  const format = $('.data_format', node).text();
  data.setFormat(asserts.assertString(format));
  if ($('input[value=\'text\']', node).is(':checked')) {
    const stringValue = $('.data_text', node).text();
    if (stringValue) {
      data.setStringValue(asserts.assertString(stringValue));
    }
  } else if ($('input[value=\'number\']', node).is(':checked')) {
    const stringValue = $('.data_text', node).text();
    const value = parseFloat(stringValue);
    if (Number.isInteger(value) && stringValue.indexOf('.') === -1) {
      data.setIntegerValue(value);
    } else if (!Number.isNaN(value)) {
      data.setFloatValue(value);
    }
  } else if ($('input[value=\'upload\']', node).is(':checked')) {
    const bytesValue = uploads.unload(node);
    if (bytesValue) {
      data.setBytesValue(bytesValue);
    }
  } else if ($('input[value=\'url\']', node).is(':checked')) {
    const url = $('.data_text', node).text();
    if (url) {
      data.setUrl(asserts.assertString(url));
    }
  }
  return data;
}
