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

goog.module('ord.identifiers');
goog.module.declareLegacyNamespace();

const asserts = goog.require('goog.asserts');

const uploads = goog.require('ord.uploads');
const utils = goog.require('ord.utils');

const ReactionIdentifier = goog.require('proto.ord.ReactionIdentifier');
const IdentifierType = goog.require('proto.ord.ReactionIdentifier.IdentifierType');

exports = {
  load,
  unload,
  add
};


/**
 * Adds and populates the reaction identifier sections in the form.
 * @param {!Array<!ReactionIdentifier>} identifiers
 */
function load(identifiers) {
  identifiers.forEach(identifier => loadIdentifier(identifier));
  if (!(identifiers.length)) {
    add();
  }
}

/**
 * Adds and populates a single reaction identifier section in the form.
 * @param {!ReactionIdentifier} identifier
 */
function loadIdentifier(identifier) {
  const node = add();
  const value = identifier.getValue();
  $('.reaction_identifier_value', node).text(value);
  utils.setSelector(node, identifier.getType());
  $('.reaction_identifier_details', node).text(identifier.getDetails());
}

/**
 * Fetches the reaction identifiers defined in the form.
 * @return {!Array<!ReactionIdentifier>}
 */
function unload() {
  const identifiers = [];
  $('.reaction_identifier').each(function(index, node) {
    node = $(node);
    if (!utils.isTemplateOrUndoBuffer(node)) {
      const identifier = unloadIdentifier(node);
      if (!utils.isEmptyMessage(identifier)) {
        identifiers.push(identifier);
      }
    }
  });
  return identifiers;
}

/**
 * Fetches a single reaction identifier defined in the form.
 * @param {!jQuery} node Root node for the identifier.
 * @return {!ReactionIdentifier}
 */
function unloadIdentifier(node) {
  const identifier = new ReactionIdentifier();

  identifier.setValue(
      asserts.assertString($('.reaction_identifier_value', node).text()));

  const type = utils.getSelectorText(node[0]);
  identifier.setType(IdentifierType[type]);
  identifier.setDetails(
      asserts.assertString($('.reaction_identifier_details', node).text()));
  return identifier;
}

/**
 * Adds a reaction identifier section to the form.
 * @return {!jQuery} The newly added parent node for the identifier.
 */
function add() {
  const node =
      utils.addSlowly('#reaction_identifier_template', $('#identifiers'));

  const uploadButton = $('.reaction_identifier_upload', node);
  uploadButton.on('change', function() {
    if ($(this).is(':checked')) {
      $('.uploader', node).show();
      $('.reaction_identifier_value', node).hide();
      $('.text_upload', node).hide();
    } else {
      $('.uploader', node).hide();
      $('.reaction_identifier_value', node).show();
    }
  });
  uploads.initialize(node);
  return node;
}
