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

goog.module('ord.outcomes');
goog.module.declareLegacyNamespace();

const asserts = goog.require('goog.asserts');

const JspbMap = goog.requireType('jspb.Map');

const data = goog.require('ord.data');
const products = goog.require('ord.products');
const utils = goog.require('ord.utils');

const Analysis = goog.require('proto.ord.Analysis');
const AnalysisType = goog.require('proto.ord.Analysis.AnalysisType');
const Data = goog.require('proto.ord.Data');
const DateTime = goog.require('proto.ord.DateTime');
const Percentage = goog.require('proto.ord.Percentage');
const ReactionOutcome = goog.require('proto.ord.ReactionOutcome');
const Time = goog.require('proto.ord.Time');

exports = {
  load,
  unload,
  add,
  addAnalysis,
  addData,
  validateOutcome,
  validateAnalysis
};


/**
 * Adds and populates the reaction outcome sections in the form.
 * @param {!Array<!ReactionOutcome>} outcomes
 */
function load(outcomes) {
  outcomes.forEach(outcome => loadOutcome(outcome));
}

/**
 * Adds and populates a reaction outcome section in the form.
 * @param {!ReactionOutcome} outcome
 */
function loadOutcome(outcome) {
  const node = add();

  const time = outcome.getReactionTime();
  if (time != null) {
    utils.writeMetric('.outcome_time', time, node);
  }
  const conversion = outcome.getConversion();
  if (conversion) {
    utils.writeMetric('.outcome_conversion', outcome.getConversion(), node);
  }

  const analysesMap = outcome.getAnalysesMap();
  analysesMap.forEach(function(analysis, name) {
    loadAnalysis(node, name, analysis);
  });

  const productsList = outcome.getProductsList();
  products.load(node, productsList);
}

/**
 * Adds and populates a reaction analysis section in the form.
 * @param {!jQuery} outcomeNode Parent reaction outcome node.
 * @param {string} name The name of this analysis.
 * @param {!Analysis} analysis
 */
function loadAnalysis(outcomeNode, name, analysis) {
  const node = addAnalysis(outcomeNode);

  $('.outcome_analysis_name', node).text(name).trigger('input');

  utils.setSelector($('.outcome_analysis_type', node), analysis.getType());
  const chmoId = analysis.getChmoId();
  if (chmoId !== 0) {
    $('.outcome_analysis_chmo_id', node).text(analysis.getChmoId());
  }
  $('.outcome_analysis_details', node).text(analysis.getDetails());

  const dataMap = analysis.getDataMap();
  dataMap.forEach(function(data, name) {
    const dataNode = addData(node);
    loadData(dataNode, name, data);
  });

  $('.outcome_analysis_manufacturer', node)
      .text(analysis.getInstrumentManufacturer());
  const calibrated = analysis.getInstrumentLastCalibrated();
  if (calibrated) {
    $('.outcome_analysis_calibrated', node).text(calibrated.getValue());
  }
  utils.setOptionalBool(
      $('.outcome_analysis_is_of_isolated_species', node),
      analysis.hasIsOfIsolatedSpecies() ? analysis.getIsOfIsolatedSpecies() :
                                          null);
}

/**
 * Adds and populates a data section in a reaction analysis.
 * @param {!jQuery} node Parent reaction analysis node.
 * @param {string} name The name of this Data record.
 * @param {!Data} dataMessage
 */
function loadData(node, name, dataMessage) {
  $('.outcome_data_name', node).text(name);
  data.loadData(node, dataMessage);
}

/**
 * Fetches the reaction outcomes defined in the form.
 * @return {!Array<!ReactionOutcome>}
 */
function unload() {
  const outcomes = [];
  $('.outcome').each(function(index, node) {
    node = $(node);
    if (!utils.isTemplateOrUndoBuffer(node)) {
      const outcome = unloadOutcome(node);
      if (!utils.isEmptyMessage(outcome)) {
        outcomes.push(outcome);
      }
    }
  });
  return outcomes;
}

/**
 * Fetches a reaction outcome defined in the form.
 * @param {!jQuery} node Root node for the reaction outcome.
 * @return {!ReactionOutcome}
 */
function unloadOutcome(node) {
  const outcome = new ReactionOutcome();

  const time = utils.readMetric('.outcome_time', new Time(), node);
  if (!utils.isEmptyMessage(time)) {
    outcome.setReactionTime(time);
  }

  const conversion =
      utils.readMetric('.outcome_conversion', new Percentage(), node);
  if (!utils.isEmptyMessage(conversion)) {
    outcome.setConversion(conversion);
  }

  const productsList = products.unload(node);
  outcome.setProductsList(productsList);

  const analysesMap = outcome.getAnalysesMap();
  $('.outcome_analysis', node).each(function(index, node) {
    node = $(node);
    if (!utils.isTemplateOrUndoBuffer(node)) {
      unloadAnalysis(node, analysesMap);
    }
  });
  return outcome;
}

/**
 * Fetches a reaction analysis defined in the form.
 * @param {!jQuery} analysisNode Root node for the reaction analysis.
 * @return {!Analysis}
 */
function unloadAnalysisSingle(analysisNode) {
  const analysis = new Analysis();
  const analysisType =
      utils.getSelectorText($('.outcome_analysis_type', analysisNode)[0]);
  analysis.setType(AnalysisType[analysisType]);
  const chmoId =
      parseInt($('.outcome_analysis_chmo_id', analysisNode).text(), 10);
  if (!isNaN(chmoId)) {
    analysis.setChmoId(chmoId);
  }
  analysis.setDetails(asserts.assertString(
      $('.outcome_analysis_details', analysisNode).text()));

  const dataMap = analysis.getDataMap();
  $('.outcome_data', analysisNode).each(function(index, dataNode) {
    dataNode = $(dataNode);
    if (!dataNode.attr('id')) {
      unloadData(dataNode, dataMap);
    }
  });
  analysis.setInstrumentManufacturer(asserts.assertString(
      $('.outcome_analysis_manufacturer', analysisNode).text()));
  const calibrated = new DateTime();
  calibrated.setValue(asserts.assertString(
      $('.outcome_analysis_calibrated', analysisNode).text()));
  if (!utils.isEmptyMessage(calibrated)) {
    analysis.setInstrumentLastCalibrated(calibrated);
  }
  const isOfIsolatedSpecies = utils.getOptionalBool(
      $('.outcome_analysis_is_of_isolated_species', analysisNode));
  if (isOfIsolatedSpecies !== null) {
    analysis.setIsOfIsolatedSpecies(isOfIsolatedSpecies);
  }

  return analysis;
}

/**
 * Fetches a reaction analysis defined in the form and adds it to `analyses`.
 * @param {!jQuery} analysisNode Root node for the reaction analysis.
 * @param {!JspbMap<string, !Analysis>} analyses
 */
function unloadAnalysis(analysisNode, analyses) {
  const analysis = unloadAnalysisSingle(analysisNode);
  const name = $('.outcome_analysis_name', analysisNode).text();
  if (name || !utils.isEmptyMessage(analysis)) {
    analyses.set(asserts.assertString(name), analysis);
  }
}

/**
 * Fetches a data record defined in the form and adds it to `dataMap`.
 * @param {!jQuery} node Root node for the Data record.
 * @param {!JspbMap<string, !Data>} dataMap
 */
function unloadData(node, dataMap) {
  const name = $('.outcome_data_name', node).text();
  const dataMessage = data.unloadData(node);
  if (name || !utils.isEmptyMessage(dataMessage)) {
    dataMap.set(asserts.assertString(name), dataMessage);
  }
}

/**
 * Adds a reaction outcome section to the form.
 * @return {!jQuery} The newly added parent node for the reaction outcome.
 */
function add() {
  const node = utils.addSlowly('#outcome_template', $('#outcomes'));
  // Add live validation handling.
  utils.addChangeHandler(node, () => { validateOutcome(node); });
  return node;
}

/**
 * Adds a reaction analysis section to the form.
 * @param {!jQuery} node Parent reaction outcome node.
 * @return {!jQuery} The newly added parent node for the reaction analysis.
 */
function addAnalysis(node) {
  const analysisNode = utils.addSlowly(
      '#outcome_analysis_template', $('.outcome_analyses', node));

  // Handle name changes.
  const nameNode = $('.outcome_analysis_name', analysisNode);
  nameNode.on('focusin', function() {
    // Store old value in val attribute.
    nameNode.data('val', nameNode.text());
  });
  nameNode.on('input', function() {
    const old_name = nameNode.data('val');
    const name = nameNode.text();
    // Remove old key.
    if (old_name) {
      // If any selector had this value selected, reset it.
      $('.analysis_key_selector', node).each(function() {
        if ($(this).val() === old_name) {
          $(this).val('');
        }
      });
      $('.analysis_key_selector option[value="' + old_name + '"]', node)
          .remove();
    }
    // Add new key.
    if (name) {
      $('.analysis_key_selector', node)
          .append('<option value="' + name + '">' + name + '</option>');
      // Ensure old value stored (necessary if focus does not change).
      nameNode.data('val', name);
    }
  });

  // Add live validation handling.
  utils.addChangeHandler(
      analysisNode, () => { validateAnalysis(analysisNode); });
  return analysisNode;
}

/**
 * Adds a new data section to the form.
 * @param {!jQuery} node Parent reaction outcome node.
 * @return {!jQuery} The newly added parent node for the Data record.
 */
function addData(node) {
  const processNode = utils.addSlowly(
      '#outcome_data_template', $('.outcome_data_repeated', node));
  data.addData(processNode);
  return processNode;
}

/**
 * Validates a reaction outcome defined in the form.
 * @param {!jQuery} node Root node for the reaction outcome.
 * @param {?jQuery=} validateNode The target node for validation results.
 */
function validateOutcome(node, validateNode = null) {
  const outcome = unloadOutcome(node);
  utils.validate(outcome, 'ReactionOutcome', node, validateNode);
}

/**
 * Validates a reaction analysis defined in the form.
 * @param {!jQuery} node Root node for the reaction analysis.
 * @param {?jQuery=} validateNode The target node for validation results.
 */
function validateAnalysis(node, validateNode = null) {
  const analysis = unloadAnalysisSingle(node);
  utils.validate(analysis, 'Analysis', node, validateNode);
}
