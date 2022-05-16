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

goog.module('ord.observations');
goog.module.declareLegacyNamespace();

const asserts = goog.require('goog.asserts');

const data = goog.require('ord.data');
const utils = goog.require('ord.utils');

const ReactionObservation = goog.require('proto.ord.ReactionObservation');
const Time = goog.require('proto.ord.Time');

exports = {
  load,
  unload,
  add,
  validateObservation
};


/**
 * Adds and populates the reaction observation sections in the form.
 * @param {!Array<!ReactionObservation>} observations
 */
function load(observations) {
  observations.forEach(observation => loadObservation(observation));
}

/**
 * Adds and populates a single reaction observation section in the form.
 * @param {!ReactionObservation} observation
 */
function loadObservation(observation) {
  const node = add();
  utils.writeMetric('.observation_time', observation.getTime(), node);
  $('.observation_comment', node).text(observation.getComment());
  data.loadData(node, observation.getImage());
}

/**
 * Fetches the reaction observations defined in the form.
 * @return {!Array<!ReactionObservation>}
 */
function unload() {
  const observations = [];
  $('.observation').each(function(index, node) {
    node = $(node);
    if (!utils.isTemplateOrUndoBuffer(node)) {
      const observation = unloadObservation(node);
      if (!utils.isEmptyMessage(observation)) {
        observations.push(observation);
      }
    }
  });
  return observations;
}

/**
 * Fetches a single reaction observation defined in the form.
 * @param {!jQuery} node Root node for the reaction observation.
 * @return {!ReactionObservation}
 */
function unloadObservation(node) {
  const observation = new ReactionObservation();
  const time = utils.readMetric('.observation_time', new Time(), node);
  if (!utils.isEmptyMessage(time)) {
    observation.setTime(time);
  }
  observation.setComment(
      asserts.assertString($('.observation_comment', node).text()));
  const image = data.unloadData(node);
  if (!utils.isEmptyMessage(image)) {
    observation.setImage(image);
  }
  return observation;
}

/**
 * Adds a reaction observation section to the form.
 * @return {!jQuery} The newly added parent node for the reaction observation.
 */
function add() {
  const node = utils.addSlowly('#observation_template', $('#observations'));
  data.addData(node);
  // Add live validation handling.
  utils.addChangeHandler(node, () => { validateObservation(node); });
  return node;
}

/**
 * Validates a single reaction observation defined in the form.
 * @param {!jQuery} node Root node for the reaction observation.
 * @param {?jQuery=} validateNode Target node for validation results.
 */
function validateObservation(node, validateNode = null) {
  const observation = unloadObservation(node);
  utils.validate(observation, 'ReactionObservation', node, validateNode);
}
