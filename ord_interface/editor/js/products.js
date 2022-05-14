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

goog.module('ord.products');
goog.module.declareLegacyNamespace();

const asserts = goog.require('goog.asserts');

const amounts = goog.require('ord.amounts');
const compounds = goog.require('ord.compounds');
const utils = goog.require('ord.utils');

const FloatValue = goog.require('proto.ord.FloatValue');
const Percentage = goog.require('proto.ord.Percentage');
const ProductCompound = goog.require('proto.ord.ProductCompound');
const Texture = goog.require('proto.ord.ProductCompound.Texture');
const TextureType = goog.require('proto.ord.ProductCompound.Texture.TextureType');
const ProductMeasurement = goog.require('proto.ord.ProductMeasurement');
const MeasurementType = goog.require('proto.ord.ProductMeasurement.MeasurementType');
const MassSpecMeasurementDetails = goog.require('proto.ord.ProductMeasurement.MassSpecMeasurementDetails');
const MassSpecMeasurementType = goog.require('proto.ord.ProductMeasurement.MassSpecMeasurementDetails.MassSpecMeasurementType');
const ReactionRoleType = goog.require('proto.ord.ReactionRole.ReactionRoleType');
const Selectivity = goog.require('proto.ord.ProductMeasurement.Selectivity');
const SelectivityType = goog.require('proto.ord.ProductMeasurement.Selectivity.SelectivityType');
const Time = goog.require('proto.ord.Time');
const Wavelength = goog.require('proto.ord.Wavelength');

exports = {
  load,
  unload,
  add,
  addMeasurement,
  validateProduct,
  validateMeasurement
};

// Freely create radio button groups by generating new input names.
let radioGroupCounter = 0;

/**
 * Adds and populates the products section of the form.
 * @param {!jQuery} node The target node for the product divs.
 * @param {!Array<!ProductCompound>} products
 */
function load(node, products) {
  products.forEach(product => loadProduct(node, product));
}

/**
 * Adds and populates a product section in the form.
 * @param {!jQuery} outcomeNode The parent ReactionOutcome node.
 * @param {!ProductCompound} product
 */
function loadProduct(outcomeNode, product) {
  const node = add(outcomeNode);

  const reactionRole = product.getReactionRole();
  utils.setSelector($('.component_reaction_role', node), reactionRole);

  const identifiers = product.getIdentifiersList();
  identifiers.forEach(identifier => {
    compounds.loadIdentifier(node, identifier);
  });

  utils.setOptionalBool(
      $('.outcome_product_desired', node),
      product.hasIsDesiredProduct() ? product.getIsDesiredProduct() : null);
  product.getMeasurementsList().forEach(
      measurement => loadMeasurement(node, measurement));
  $('.outcome_product_color', node).text(product.getIsolatedColor());
  const texture = product.getTexture();
  if (texture) {
    utils.setSelector(
        $('.outcome_product_texture_type', node), texture.getType());
    $('.outcome_product_texture_details', node).text(texture.getDetails());
  }
  const features = product.getFeaturesMap();
  features.forEach(function(feature, name) {
    const featureNode = compounds.addFeature(node);
    compounds.loadFeature(featureNode, name, feature);
  });
}

/**
 * Fetches the products defined in the form.
 * @param {!jQuery} node The parent ReactionOutcome node.
 * @return {!Array<!ProductCompound>}
 */
function unload(node) {
  const products = [];
  $('.outcome_product', node).each(function(index, productNode) {
    productNode = $(productNode);
    if (!productNode.attr('id')) {
      // Not a template.
      const product = unloadProduct(productNode);
      if (!utils.isEmptyMessage(product)) {
        products.push(product);
      }
    }
  });
  return products;
}

/**
 * Fetches a product defined in the form.
 * @param {!jQuery} node An element containing a product.
 * @return {!ProductCompound}
 */
function unloadProduct(node) {
  const product = new ProductCompound();

  const reactionRole =
      utils.getSelectorText($('.component_reaction_role', node)[0]);
  product.setReactionRole(ReactionRoleType[reactionRole]);

  const identifiers =
      compounds.unloadIdentifiers($('.product_compound_identifiers', node));

  if (identifiers.some(e => !utils.isEmptyMessage(e))) {
    product.setIdentifiersList(identifiers);
  }

  const isDesiredProduct =
      utils.getOptionalBool($('.outcome_product_desired', node));
  if (isDesiredProduct !== null) {
    product.setIsDesiredProduct(isDesiredProduct);
  }

  const measurements = [];
  $('.product_measurement', node).each(function(index, measurementNode) {
    measurementNode = $(measurementNode);
    if (!measurementNode.attr('id')) {
      // Not a template.
      const measurement = unloadMeasurement(measurementNode);
      if (!utils.isEmptyMessage(measurement)) {
        measurements.push(measurement);
      }
    }
  });
  product.setMeasurementsList(measurements);

  product.setIsolatedColor(
      asserts.assertString($('.outcome_product_color', node).text()));

  const texture = new Texture();
  const textureType =
      utils.getSelectorText($('.outcome_product_texture_type', node)[0]);
  texture.setType(TextureType[textureType]);
  texture.setDetails(
      asserts.assertString($('.outcome_product_texture_details', node).text()));
  if (!utils.isEmptyMessage(texture)) {
    product.setTexture(texture);
  }

  const featuresMap = product.getFeaturesMap();
  $('.feature', node).each(function(index, featureNode) {
    featureNode = $(featureNode);
    if (!featureNode.attr('id')) {
      compounds.unloadFeature(featureNode, featuresMap);
    }
  });

  return product;
}

/**
 * Fetches the set of Analysis keys.
 * @param {!jQuery} node The parent ReactionOutcome div.
 * @param {string} tag Analysis target, e.g. "identities", "yields", etc.
 * @return {!Array<string>}
 */
function unloadAnalysisKeys(node, tag) {
  const values = [];
  $('.outcome_product_analysis_' + tag, node).each(function(index, tagNode) {
    tagNode = $(tagNode);
    if (!tagNode.attr('id')) {
      // Not a template.
      const value = $('.analysis_key_selector', tagNode).val();
      if (value !== '') {
        values.push(value);
      }
    }
  });
  return values;
}

/**
 * Adds a reaction product section to the form.
 * @param {!jQuery} node Target ReactionOutcome node for the new product.
 * @return {!jQuery} The newly created node.
 */
function add(node) {
  const productNode = utils.addSlowly(
      '#outcome_product_template', $('.outcome_products', node));
  amounts.init(node);

  // Connect reaction role selection to is_desired_product.
  const roleSelector = $('.component_reaction_role', node);
  roleSelector.on('change', function() {
    if (utils.getSelectorText(this) === 'PRODUCT') {
      $('.is_desired_product', node).show();
    } else {
      $('.is_desired_product', node).hide();
    }
  });
  roleSelector.trigger('change');

  // Add live validation handling.
  utils.addChangeHandler(productNode, () => {
    validateProduct(productNode);
  });
  return productNode;
}

/**
 * Adds keys for defined analyses to the analysis selector.
 * @param {!jQuery} node Parent node containing ReactionOutcome data.
 * @param {!jQuery} analysisSelectorNode Node containing an analysis selector.
 */
function populateAnalysisSelector(node, analysisSelectorNode) {
  const outcomeNode = node.closest('.outcome');
  $('.outcome_analysis_name', outcomeNode).each(function() {
    const name = $(this).text();
    if (name) {
      $('.analysis_key_selector', analysisSelectorNode)
          .append('<option value="' + name + '">' + name + '</option>');
    }
  });
}

/**
 * Adds a ProductMeasurement section to the form.
 * @param {!jQuery} node Parent node for the ProductCompound.
 * @return {!jQuery} The newly created node.
 */
function addMeasurement(node) {
  const measurementNode = utils.addSlowly(
      '#product_measurement_template',
      $('.product_measurement_repeated', node));
  amounts.init(measurementNode);
  populateAnalysisSelector(
      node, $('.product_measurement_analysis_key', measurementNode));

  // Set up the radio buttons for the value type.
  const buttons = $('.product_measurement_value_type input', measurementNode);
  buttons.attr('name', 'product_measurements_' + radioGroupCounter++);
  buttons.on('change', function() {
    if (this.value === 'string') {
      $('.product_measurement_pm', measurementNode).hide();
      $('.product_measurement_precision', measurementNode).hide();
    } else {
      $('.product_measurement_pm', measurementNode).show();
      $('.product_measurement_precision', measurementNode).show();
    }
    if (this.value === 'mass') {
      $('.amount_units', measurementNode).show();
    } else {
      $('.amount_units', measurementNode).hide();
    }
  });

  // Add an empty compound node for the authentic standard.
  const authenticStandard = compounds.add(measurementNode);
  const title = $('.h4', authenticStandard);
  title.text('Authentic Standard');
  // Remove the "amount" section from the authentic standard.
  $('.amount', authenticStandard).hide();
  // Adjust the heading sizes.
  $('.h4', authenticStandard).addClass('h7').removeClass('h4');
  $('.h5', authenticStandard).addClass('h8').removeClass('h5');

  // Show/hide the authentic standard based on the optional bool.
  const usesAuthenticStandard =
      $('.product_measurement_uses_authentic_standard', measurementNode);
  usesAuthenticStandard.on('change', function() {
    if (utils.getOptionalBool(usesAuthenticStandard)) {
      $('.product_measurement_authentic_standard', measurementNode).show();
    } else {
      $('.product_measurement_authentic_standard', measurementNode).hide();
    }
  });
  usesAuthenticStandard.trigger('change');

  // Add live validation handling.
  utils.addChangeHandler(measurementNode, () => {
    validateMeasurement(measurementNode);
  });

  // Show/hide fields based on the measurement type.
  const measurementTypeSelector =
      $('.product_measurement_type select', measurementNode);
  measurementTypeSelector.on('change', function() {
    const measurementType = this.options[this.selectedIndex].text;
    if (measurementType === 'UNSPECIFIED') {
      $('.product_measurement_value_group', measurementNode).hide();
      $('.retention_time', measurementNode).hide();
      $('.wavelength', measurementNode).hide();
      $('.mass_spec_details', measurementNode).hide();
      $('.selectivity', measurementNode).hide();
    } else if (measurementType === 'CUSTOM') {
      $('.product_measurement_value_group', measurementNode).show();
      $('.retention_time', measurementNode).show();
      $('.wavelength', measurementNode).show();
      $('.mass_spec_details', measurementNode).show();
      $('.selectivity', measurementNode).show();
    } else if (measurementType === 'IDENTITY') {
      $('.product_measurement_value_group', measurementNode).hide();
      $('.retention_time', measurementNode).show();
      $('.wavelength', measurementNode).hide();
      $('.mass_spec_details', measurementNode).hide();
      $('.selectivity', measurementNode).hide();
    } else if (measurementType === 'YIELD') {
      $('.product_measurement_value_group', measurementNode).show();
      $('.retention_time', measurementNode).hide();
      $('.wavelength', measurementNode).hide();
      $('.mass_spec_details', measurementNode).hide();
      $('.selectivity', measurementNode).hide();
      $('.product_measurement_percentage', measurementNode).trigger('click');
    } else if (measurementType === 'SELECTIVITY') {
      $('.product_measurement_value_group', measurementNode).show();
      $('.retention_time', measurementNode).hide();
      $('.wavelength', measurementNode).hide();
      $('.mass_spec_details', measurementNode).hide();
      $('.selectivity', measurementNode).show();
      $('.product_measurement_string', measurementNode).trigger('click');
    } else if (measurementType === 'PURITY') {
      $('.product_measurement_value_group', measurementNode).show();
      $('.retention_time', measurementNode).hide();
      $('.wavelength', measurementNode).show();
      $('.mass_spec_details', measurementNode).hide();
      $('.selectivity', measurementNode).hide();
      $('.product_measurement_percentage', measurementNode).trigger('click');
    } else if (measurementType === 'AREA') {
      $('.product_measurement_value_group', measurementNode).show();
      $('.retention_time', measurementNode).show();
      $('.wavelength', measurementNode).show();
      $('.mass_spec_details', measurementNode).show();
      $('.selectivity', measurementNode).hide();
      $('.product_measurement_float', measurementNode).trigger('click');
    } else if (measurementType === 'COUNTS') {
      $('.product_measurement_value_group', measurementNode).show();
      $('.retention_time', measurementNode).show();
      $('.wavelength', measurementNode).hide();
      $('.mass_spec_details', measurementNode).show();
      $('.selectivity', measurementNode).hide();
      $('.product_measurement_float', measurementNode).trigger('click');
    } else if (measurementType === 'INTENSITY') {
      $('.product_measurement_value_group', measurementNode).show();
      $('.retention_time', measurementNode).show();
      $('.wavelength', measurementNode).show();
      $('.mass_spec_details', measurementNode).show();
      $('.selectivity', measurementNode).hide();
      $('.product_measurement_float', measurementNode).trigger('click');
    } else if (measurementType === 'AMOUNT') {
      $('.product_measurement_value_group', measurementNode).show();
      $('.retention_time', measurementNode).hide();
      $('.wavelength', measurementNode).hide();
      $('.mass_spec_details', measurementNode).hide();
      $('.selectivity', measurementNode).hide();
      $('.product_measurement_mass', measurementNode).trigger('click');
    }
  });
  measurementTypeSelector.trigger('change');

  return measurementNode;
}

/**
 * Adds and populates a ProductMeasurement section in the form.
 * @param {!jQuery} productNode The parent ProductCompound node.
 * @param {!ProductMeasurement} measurement
 */
function loadMeasurement(productNode, measurement) {
  const node = addMeasurement(productNode);
  $('.analysis_key_selector', node).val(measurement.getAnalysisKey());
  utils.setSelector(
      $('.product_measurement_type', node), measurement.getType());
  $('.product_measurement_type select', node).trigger('change');
  $('.product_measurement_details', node).text(measurement.getDetails());
  utils.setOptionalBool(
      $('.product_measurement_uses_internal_standard', node),
      measurement.hasUsesInternalStandard() ?
          measurement.getUsesInternalStandard() :
          null);
  utils.setOptionalBool(
      $('.product_measurement_is_normalized', node),
      measurement.hasIsNormalized() ? measurement.getIsNormalized() : null);
  utils.setOptionalBool(
      $('.product_measurement_uses_authentic_standard', node),
      measurement.hasUsesAuthenticStandard() ?
          measurement.getUsesAuthenticStandard() :
          null);

  const authenticStandard = measurement.getAuthenticStandard();
  if (authenticStandard) {
    compounds.loadIntoCompound(node, authenticStandard);
  }
  $('.product_measurement_uses_authentic_standard', node).trigger('change');

  if (measurement.hasPercentage()) {
    $('.product_measurement_percentage', node).trigger('click');
    $('.product_measurement_value', node).addClass('floattext');
    if (measurement.getPercentage().hasValue()) {
      $('.product_measurement_value', node)
          .text(measurement.getPercentage().getValue());
    }
    if (measurement.getPercentage().hasPrecision()) {
      $('.product_measurement_precision', node)
          .text(measurement.getPercentage().getPrecision());
    }
  } else if (measurement.hasFloatValue()) {
    $('.product_measurement_float', node).trigger('click');
    $('.product_measurement_value', node).addClass('floattext');
    if (measurement.getFloatValue().hasValue()) {
      $('.product_measurement_value', node)
          .text(measurement.getFloatValue().getValue());
    }
    if (measurement.getFloatValue().hasPrecision()) {
      $('.product_measurement_precision', node)
          .text(measurement.getFloatValue().getPrecision());
    }
  } else if (measurement.getStringValue()) {
    $('.product_measurement_string', node).trigger('click');
    $('.product_measurement_value', node).removeClass('floattext');
    if (measurement.getStringValue()) {
      $('.product_measurement_value', node).text(measurement.getStringValue());
    }
  } else if (measurement.hasAmount()) {
    $('.product_measurement_mass', node).trigger('click');
    const valueNode = $('.product_measurement_value_type', node);
    amounts.load(valueNode, measurement.getAmount());
  }

  const retentionTime = measurement.getRetentionTime();
  if (retentionTime) {
    utils.writeMetric(
        '.product_measurement_retention_time', retentionTime, node);
  }

  const massSpec = measurement.getMassSpecDetails();
  if (massSpec) {
    utils.setSelector(
        $('.product_measurement_mass_spec_type', node), massSpec.getType());
    $('.product_measurement_mass_spec_details', node)
        .text(massSpec.getDetails());
    if (massSpec.hasTicMinimumMz()) {
      $('.product_measurement_mass_spec_tic_minimum_mz', node)
          .text(massSpec.getTicMinimumMz());
    }
    if (massSpec.hasTicMaximumMz()) {
      $('.product_measurement_mass_spec_tic_maximum_mz', node)
          .text(massSpec.getTicMaximumMz());
    }
    const eicMasses = massSpec.getEicMassesList().join(',');
    $('.product_measurement_mass_spec_eic_masses', node).text(eicMasses);
  }

  const selectivity = measurement.getSelectivity();
  if (selectivity) {
    utils.setSelector(
        $('.product_measurement_selectivity_type', node),
        selectivity.getType());
    $('.product_measurement_selectivity_details', node)
        .text(selectivity.getDetails());
  }

  const wavelength = measurement.getWavelength();
  if (wavelength) {
    utils.writeMetric('.product_measurement_wavelength', wavelength, node);
  }
}

/**
 * Fetches a ProductMeasurement defined in the form.
 * @param {!jQuery} node An element containing a ProductMeasurement.
 * @return {!ProductMeasurement}
 */
function unloadMeasurement(node) {
  const measurement = new ProductMeasurement();
  const analysisKey = $('.product_measurement_analysis_key select', node).val();
  if (analysisKey) {
    measurement.setAnalysisKey(asserts.assertString(analysisKey));
  }
  const measurementType =
      utils.getSelectorText($('.product_measurement_type', node)[0]);
  measurement.setType(MeasurementType[measurementType]);
  measurement.setDetails(
      asserts.assertString($('.product_measurement_details', node).text()));
  const usesInternalStandard = utils.getOptionalBool(
      $('.product_measurement_uses_internal_standard', node));
  if (usesInternalStandard !== null) {
    measurement.setUsesInternalStandard(usesInternalStandard);
  }
  const isNormalized =
      utils.getOptionalBool($('.product_measurement_is_normalized', node));
  if (isNormalized !== null) {
    measurement.setIsNormalized(isNormalized);
  }
  const usesAuthenticStandard = utils.getOptionalBool(
      $('.product_measurement_uses_authentic_standard', node));
  if (usesAuthenticStandard !== null) {
    measurement.setUsesAuthenticStandard(usesAuthenticStandard);
  }

  const authenticStandardNode =
      $('.product_measurement_authentic_standard', node);
  const compound = compounds.unloadCompound(authenticStandardNode);
  if (!utils.isEmptyMessage(compound)) {
    measurement.setAuthenticStandard(compound);
  }

  if ($('.product_measurement_percentage', node).is(':checked')) {
    const percentage = new Percentage();
    const value = parseFloat($('.product_measurement_value', node).text());
    if (!isNaN(value)) {
      percentage.setValue(value);
    }
    const precision =
        parseFloat($('.product_measurement_precision', node).text());
    if (!isNaN(precision)) {
      percentage.setPrecision(precision);
    }
    if (!utils.isEmptyMessage(percentage)) {
      measurement.setPercentage(percentage);
    }
  } else if ($('.product_measurement_float', node).is(':checked')) {
    const floatValue = new FloatValue();
    const value = parseFloat($('.product_measurement_value', node).text());
    if (!isNaN(value)) {
      floatValue.setValue(value);
    }
    const precision =
        parseFloat($('.product_measurement_precision', node).text());
    if (!isNaN(precision)) {
      floatValue.setPrecision(precision);
    }
    if (!utils.isEmptyMessage(floatValue)) {
      measurement.setFloatValue(floatValue);
    }
  } else if ($('.product_measurement_string', node).is(':checked')) {
    const stringValue = $('.product_measurement_value', node).text();
    if (stringValue) {
      measurement.setStringValue(asserts.assertString(stringValue));
    }
  } else if ($('.product_measurement_mass', node).is(':checked')) {
    const amount = amounts.unload($('.product_measurement_value_type', node));
    if (!utils.isEmptyMessage(amount)) {
      measurement.setAmount(amount);
    }
  }

  const retentionTime =
      utils.readMetric('.product_measurement_retention_time', new Time(), node);
  if (!utils.isEmptyMessage(retentionTime)) {
    measurement.setRetentionTime(retentionTime);
  }

  const massSpecDetails = new MassSpecMeasurementDetails();
  const massSpecType =
      utils.getSelectorText($('.product_measurement_mass_spec_type', node)[0]);
  massSpecDetails.setType(MassSpecMeasurementType[massSpecType]);
  massSpecDetails.setDetails(asserts.assertString(
      $('.product_measurement_mass_spec_details', node).text()));
  const ticMinimumMz = parseFloat(
      $('.product_measurement_mass_spec_tic_minimum_mz', node).text());
  if (!isNaN(ticMinimumMz)) {
    massSpecDetails.setTicMinimumMz(ticMinimumMz);
  }
  const ticMaximumMz = parseFloat(
      $('.product_measurement_mass_spec_tic_maximum_mz', node).text());
  if (!isNaN(ticMaximumMz)) {
    massSpecDetails.setTicMaximumMz(ticMaximumMz);
  }
  if ($('.product_measurement_mass_spec_eic_masses', node).text()) {
    const eicMasses = $('.product_measurement_mass_spec_eic_masses', node)
                          .text()
                          .split(',')
                          .map(parseFloat);
    massSpecDetails.setEicMassesList(eicMasses);
  }
  if (!utils.isEmptyMessage(massSpecDetails)) {
    measurement.setMassSpecDetails(massSpecDetails);
  }

  const selectivity = new Selectivity();
  const selectivityType = utils.getSelectorText(
      $('.product_measurement_selectivity_type', node)[0]);
  selectivity.setType(SelectivityType[selectivityType]);
  selectivity.setDetails(asserts.assertString(
      $('.product_measurement_selectivity_details', node).text()));
  if (!utils.isEmptyMessage(selectivity)) {
    measurement.setSelectivity(selectivity);
  }

  const wavelength = utils.readMetric(
      '.product_measurement_wavelength', new Wavelength(), node);
  if (!utils.isEmptyMessage(wavelength)) {
    measurement.setWavelength(wavelength);
  }
  return measurement;
}

/**
 * Validates a ProductMeasurement defined in the form.
 * @param {!jQuery} node A node containing a ProductMeasurement.
 * @param {?jQuery=} validateNode The target div for validation results.
 */
function validateMeasurement(node, validateNode = null) {
  const measurement = unloadMeasurement(node);
  utils.validate(measurement, 'ProductMeasurement', node, validateNode);
}

/**
 * Validates a product defined in the form.
 * @param {!jQuery} node A node containing a reaction product.
 * @param {?jQuery=} validateNode The target div for validation results.
 */
function validateProduct(node, validateNode = null) {
  const product = unloadProduct(node);
  utils.validate(product, 'ProductCompound', node, validateNode);
  compounds.renderCompound(node, product);
}
