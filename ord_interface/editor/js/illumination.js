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

goog.module('ord.illumination');
goog.module.declareLegacyNamespace();

const asserts = goog.require('goog.asserts');

const utils = goog.require('ord.utils');

const IlluminationConditions = goog.require('proto.ord.IlluminationConditions');
const IlluminationType = goog.require('proto.ord.IlluminationConditions.IlluminationType');
const Length = goog.require('proto.ord.Length');
const Wavelength = goog.require('proto.ord.Wavelength');

exports = {
  load,
  unload,
  validateIllumination
};


/**
 * Adds and populates the illumination conditions section in the form.
 * @param {!IlluminationConditions} illumination
 */
function load(illumination) {
  utils.setSelector($('#illumination_type'), illumination.getType());
  $('#illumination_details').text(illumination.getDetails());
  const wavelength = illumination.getPeakWavelength();
  utils.writeMetric('#illumination_wavelength', wavelength);
  $('#illumination_color').text(illumination.getColor());
  const distance = illumination.getDistanceToVessel();
  utils.writeMetric('#illumination_distance', distance);
}

/**
 * Fetches the illumination conditions defined in the form.
 * @return {!IlluminationConditions}
 */
function unload() {
  const illumination = new IlluminationConditions();
  const illuminationType = utils.getSelectorText($('#illumination_type')[0]);
  illumination.setType(IlluminationType[illuminationType]);
  illumination.setDetails(
      asserts.assertString($('#illumination_details').text()));

  const wavelength =
      utils.readMetric('#illumination_wavelength', new Wavelength());
  if (!utils.isEmptyMessage(wavelength)) {
    illumination.setPeakWavelength(wavelength);
  }
  illumination.setColor(asserts.assertString($('#illumination_color').text()));
  const distance = utils.readMetric('#illumination_distance', new Length());
  if (!utils.isEmptyMessage(distance)) {
    illumination.setDistanceToVessel(distance);
  }
  return illumination;
}

/**
 * Validates the illumination conditions defined in the form.
 * @param {!jQuery} node Root node for the illumination conditions.
 * @param {?jQuery=} validateNode Target node for validation results.
 */
function validateIllumination(node, validateNode = null) {
  const illumination = unload();
  utils.validate(illumination, 'IlluminationConditions', node, validateNode);
}
