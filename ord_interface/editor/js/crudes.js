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

goog.module('ord.crudes');
goog.module.declareLegacyNamespace();
exports = {
  load,
  unload,
  add
};

const asserts = goog.require('goog.asserts');

const amounts = goog.require('ord.amounts');
const utils = goog.require('ord.utils');

const CrudeComponent = goog.require('proto.ord.CrudeComponent');

/**
 * Adds and populates the crude components of a reaction input.
 * @param {!jQuery} node Root node for the parent reaction input.
 * @param {!Array<!CrudeComponent>} crudes
 */
function load(node, crudes) {
  crudes.forEach(crude => loadCrude(node, crude));
}

/**
 * Adds and populates a single crude component section in the form.
 * @param {!jQuery} root Root node for the parent reaction input.
 * @param {!CrudeComponent} crude
 */
function loadCrude(root, crude) {
  const node = add(root);

  const reactionId = crude.getReactionId();
  $('.crude_reaction', node).text(reactionId);

  const workup = crude.hasIncludesWorkup() ? crude.getIncludesWorkup() : null;
  utils.setOptionalBool($('.crude_includes_workup', node), workup);

  const derived =
      crude.hasHasDerivedAmount() ? crude.getHasDerivedAmount() : null;
  utils.setOptionalBool($('.crude_has_derived', node), derived);

  const amount = crude.getAmount();
  amounts.load(node, amount);
}

/**
 * Fetches the crude components defined for a reaction input in the form.
 * @param {!jQuery} node Root node for the parent reaction input.
 * @return {!Array<!CrudeComponent>}
 */
function unload(node) {
  const crudes = [];
  $('.crude', node).each(function(index, crudeNode) {
    crudeNode = $(crudeNode);
    if (!crudeNode.attr('id')) {
      // Not a template.
      const crude = unloadCrude(crudeNode);
      if (!utils.isEmptyMessage(crude)) {
        crudes.push(crude);
      }
    }
  });
  return crudes;
}

/**
 * Fetches a single crude component defined in the form.
 * @param {!jQuery} node Root node for the crude component.
 * @return {!CrudeComponent}
 */
function unloadCrude(node) {
  const crude = new CrudeComponent();

  const reactionId = $('.crude_reaction', node).text();
  crude.setReactionId(asserts.assertString(reactionId));

  const workup = utils.getOptionalBool($('.crude_includes_workup', node));
  if (workup !== null) {
    crude.setIncludesWorkup(workup);
  }

  const derived = utils.getOptionalBool($('.crude_has_derived', node));
  if (derived !== null) {
    crude.setHasDerivedAmount(derived);
  }

  const amount = amounts.unload(node);
  if (!utils.isEmptyMessage(amount)) {
    crude.setAmount(amount);
  }

  return crude;
}

/**
 * Adds a crude component section to the given reaction input.
 * @param {!jQuery} root Root node for the parent reaction input.
 * @return {!jQuery} The newly added root node for the crude component.
 */
function add(root) {
  const node = utils.addSlowly('#crude_template', $('.crudes', root));
  amounts.init(node);
  return node;
}
