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

goog.module('ord.notes');
goog.module.declareLegacyNamespace();

const asserts = goog.require('goog.asserts');

const utils = goog.require('ord.utils');

const ReactionNotes = goog.require('proto.ord.ReactionNotes');

exports = {
  load,
  unload,
  validateNotes
};


/**
 * Adds and populates the reaction nodes section in the form.
 * @param {!ReactionNotes} notes
 */
function load(notes) {
  utils.setOptionalBool(
      $('#notes_heterogeneous'),
      notes.hasIsHeterogeneous() ? notes.getIsHeterogeneous() : null);
  utils.setOptionalBool(
      $('#notes_precipitate'),
      notes.hasFormsPrecipitate() ? notes.getFormsPrecipitate() : null);
  utils.setOptionalBool(
      $('#notes_exothermic'),
      notes.hasIsExothermic() ? notes.getIsExothermic() : null);
  utils.setOptionalBool(
      $('#notes_offgas'), notes.hasOffgasses() ? notes.getOffgasses() : null);
  utils.setOptionalBool(
      $('#notes_moisture'),
      notes.hasIsSensitiveToMoisture() ? notes.getIsSensitiveToMoisture() :
                                         null);
  utils.setOptionalBool(
      $('#notes_oxygen'),
      notes.hasIsSensitiveToOxygen() ? notes.getIsSensitiveToOxygen() : null);
  utils.setOptionalBool(
      $('#notes_light'),
      notes.hasIsSensitiveToLight() ? notes.getIsSensitiveToLight() : null);
  $('#notes_safety').text(notes.getSafetyNotes());
  $('#notes_details').text(notes.getProcedureDetails());
}

/**
 * Fetches the reaction notes defined in the form.
 * @return {!ReactionNotes}
 */
function unload() {
  const notes = new ReactionNotes();
  const isHeterogeneous = utils.getOptionalBool($('#notes_heterogeneous'));
  if (isHeterogeneous !== null) {
    notes.setIsHeterogeneous(isHeterogeneous);
  }
  const formsPrecipitate = utils.getOptionalBool($('#notes_precipitate'));
  if (formsPrecipitate !== null) {
    notes.setFormsPrecipitate(formsPrecipitate);
  }
  const isExothermic = utils.getOptionalBool($('#notes_exothermic'));
  if (isExothermic !== null) {
    notes.setIsExothermic(isExothermic);
  }
  const offgasses = utils.getOptionalBool($('#notes_offgas'));
  if (offgasses !== null) {
    notes.setOffgasses(offgasses);
  }
  const isSensitiveToMoisture = utils.getOptionalBool($('#notes_moisture'));
  if (isSensitiveToMoisture !== null) {
    notes.setIsSensitiveToMoisture(isSensitiveToMoisture);
  }
  const isSensitiveToOxygen = utils.getOptionalBool($('#notes_oxygen'));
  if (isSensitiveToOxygen !== null) {
    notes.setIsSensitiveToOxygen(isSensitiveToOxygen);
  }
  const isSensitiveToLight = utils.getOptionalBool($('#notes_light'));
  if (isSensitiveToLight !== null) {
    notes.setIsSensitiveToLight(isSensitiveToLight);
  }
  notes.setSafetyNotes(asserts.assertString($('#notes_safety').text()));
  notes.setProcedureDetails(asserts.assertString($('#notes_details').text()));
  return notes;
}

/**
 * Validates the reaction notes defined in the form.
 * @param {!jQuery} node Root node for the reaction notes.
 * @param {?jQuery=} validateNode Target node for validation results.
 */
function validateNotes(node, validateNode = null) {
  const notes = unload();
  utils.validate(notes, 'ReactionNotes', node, validateNode);
}
