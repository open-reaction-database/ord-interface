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

goog.module('ord.stirring');
goog.module.declareLegacyNamespace();

const asserts = goog.require('goog.asserts');

const utils = goog.require('ord.utils');

const StirringConditions = goog.require('proto.ord.StirringConditions');
const StirringMethodType = goog.require('proto.ord.StirringConditions.StirringMethodType');
const StirringRate = goog.require('proto.ord.StirringConditions.StirringRate');
const StirringRateType = goog.require('proto.ord.StirringConditions.StirringRate.StirringRateType');

exports = {
  load,
  unload,
  validateStirring
};


/**
 * Adds and populates the stirring conditions section in the form.
 * @param {!StirringConditions} stirring
 * @param {!jQuery} target
 */
function load(stirring, target) {
  utils.setSelector($('.stirring_method_type', target), stirring.getType());
  $('.stirring_method_details', target).text(stirring.getDetails());
  const rate = stirring.getRate();
  if (rate) {
    utils.setSelector($('.stirring_rate_type', target), rate.getType());
    $('.stirring_rate_details', target).text(rate.getDetails());
    const rpm = rate.getRpm();
    if (rpm !== 0) {
      $('.stirring_rpm', target).text(rpm);
    }
  }
}

/**
 * Fetches the stirring conditions from the form.
 * @param {!jQuery} target
 * @return {!StirringConditions}
 */
function unload(target) {
  const stirring = new StirringConditions();
  const methodType =
      utils.getSelectorText($('.stirring_method_type', target)[0]);
  stirring.setType(StirringMethodType[methodType]);
  stirring.setDetails(
      asserts.assertString($('.stirring_method_details', target).text()));

  const rate = new StirringRate();
  const rateType = utils.getSelectorText($('.stirring_rate_type', target)[0]);
  rate.setType(StirringRateType[rateType]);
  rate.setDetails(
      asserts.assertString($('.stirring_rate_details', target).text()));
  const rpm = parseFloat($('.stirring_rpm', target).text());
  if (!isNaN(rpm)) {
    rate.setRpm(rpm);
  }
  if (!utils.isEmptyMessage(rate)) {
    stirring.setRate(rate);
  }
  return stirring;
}

/**
 * Validates the stirring conditions defined in the form.
 * @param {!jQuery} node The div containing the stirring conditions.
 * @param {?jQuery=} validateNode The target div for validation results.
 */
function validateStirring(node, validateNode = null) {
  const stirring = unload(node);
  utils.validate(stirring, 'StirringConditions', node, validateNode);
}
