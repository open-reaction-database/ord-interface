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

goog.module('ord.provenance');
goog.module.declareLegacyNamespace();

const asserts = goog.require('goog.asserts');

const utils = goog.require('ord.utils');

const DateTime = goog.require('proto.ord.DateTime');
const Person = goog.require('proto.ord.Person');
const ReactionProvenance = goog.require('proto.ord.ReactionProvenance');
const RecordEvent = goog.require('proto.ord.RecordEvent');

exports = {
  load,
  unload,
  addModification,
  validateProvenance
};


/**
 * Adds and populates the provenance section in the form.
 * @param {!ReactionProvenance} provenance
 */
function load(provenance) {
  const experimenter = provenance.getExperimenter();
  if (experimenter) {
    loadPerson($('#provenance_experimenter'), experimenter);
  }
  $('#provenance_city').text(provenance.getCity());

  const start = provenance.getExperimentStart();
  if (start) {
    $('#provenance_start').text(start.getValue());
  }

  $('#provenance_doi').text(provenance.getDoi());
  $('#provenance_patent').text(provenance.getPatent());
  $('#provenance_url').text(provenance.getPublicationUrl());
  loadRecordEvent($('#provenance_created'), provenance.getRecordCreated());
  provenance.getRecordModifiedList().forEach(modified => {
    const node = addModification();
    loadRecordEvent(node, modified);
  });
}

/**
 * Adds and populates a record event section in the form.
 * @param {!jQuery} node The target div.
 * @param {?RecordEvent} record
 */
function loadRecordEvent(node, record) {
  if (!record) {
    return;
  }
  const time = record.getTime();
  if (time) {
    $('.provenance_time', node).text(time.getValue());
  }
  loadPerson(node, record.getPerson());
  $('.provenance_details', node).text(record.getDetails());
}

/**
 * Adds and populates a person section in the form.
 * @param {!jQuery} node The target div.
 * @param {?Person} person
 */
function loadPerson(node, person) {
  if (!person) {
    return;
  }
  $('.provenance_username', node).text(person.getUsername());
  $('.provenance_name', node).text(person.getName());
  $('.provenance_orcid', node).text(person.getOrcid());
  $('.provenance_organization', node).text(person.getOrganization());
  $('.provenance_email', node).text(person.getEmail());
}

/**
 * Fetches reaction provenance as defined in the form.
 * @return {!ReactionProvenance}
 */
function unload() {
  const provenance = new ReactionProvenance();

  const experimenter = unloadPerson($('#provenance_experimenter'));
  if (!utils.isEmptyMessage(experimenter)) {
    provenance.setExperimenter(experimenter);
  }

  provenance.setCity(asserts.assertString($('#provenance_city').text()));

  const start = new DateTime();
  start.setValue(asserts.assertString($('#provenance_start').text()));
  if (!utils.isEmptyMessage(start)) {
    provenance.setExperimentStart(start);
  }

  provenance.setDoi(asserts.assertString($('#provenance_doi').text()));
  provenance.setPatent(asserts.assertString($('#provenance_patent').text()));
  provenance.setPublicationUrl(
      asserts.assertString($('#provenance_url').text()));

  const created = unloadRecordEvent($('#provenance_created'));
  if (!utils.isEmptyMessage(created)) {
    provenance.setRecordCreated(created);
  }

  const modifieds = [];
  $('.provenance_modified', $('#provenance_modifieds'))
      .each(function(index, node) {
        node = $(node);
        if (!utils.isTemplateOrUndoBuffer(node)) {
          const modified = unloadRecordEvent(node);
          if (!utils.isEmptyMessage(modified)) {
            modifieds.push(modified);
          }
        }
      });
  provenance.setRecordModifiedList(modifieds);
  return provenance;
}

/**
 * Fetches a record event as defined in the form.
 * @param {!jQuery} node Parent node containing the record event.
 * @return {!RecordEvent}
 */
function unloadRecordEvent(node) {
  const created = new RecordEvent();
  const createdTime = new DateTime();
  createdTime.setValue(
      asserts.assertString($('.provenance_time', node).text()));
  if (!utils.isEmptyMessage(createdTime)) {
    created.setTime(createdTime);
  }
  const createdPerson = unloadPerson(node);
  if (!utils.isEmptyMessage(createdPerson)) {
    created.setPerson(createdPerson);
  }
  created.setDetails(
      asserts.assertString($('.provenance_details', node).text()));
  return created;
}

/**
 * Fetches a person message as defined in the form.
 * @param {!jQuery} node Parent node containing the person message.
 * @return {!Person}
 */
function unloadPerson(node) {
  const person = new Person();
  person.setUsername(
      asserts.assertString($('.provenance_username', node).text()));
  person.setName(asserts.assertString($('.provenance_name', node).text()));
  person.setOrcid(asserts.assertString($('.provenance_orcid', node).text()));
  person.setOrganization(
      asserts.assertString($('.provenance_organization', node).text()));
  person.setEmail(asserts.assertString($('.provenance_email', node).text()));
  return person;
}

/**
 * Adds a record_modified section to the form.
 * @return {!jQuery} The div containing the new event.
 */
function addModification() {
  return utils.addSlowly(
      '#provenance_modified_template', $('#provenance_modifieds'));
}

/**
 * Validates the reaction provenence defined in the form.
 * @param {!jQuery} node The node containing reaction provenance information.
 * @param {?jQuery=} validateNode The target div for validation results.
 */
function validateProvenance(node, validateNode = null) {
  const provenance = unload();
  utils.validate(provenance, 'ReactionProvenance', node, validateNode);
}
