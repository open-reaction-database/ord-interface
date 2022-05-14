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

goog.module('ord.amounts');
goog.module.declareLegacyNamespace();
exports = {
  init,
  load,
  unload,
};

const utils = goog.require('ord.utils');

const Amount = goog.require('proto.ord.Amount');
const Mass = goog.require('proto.ord.Mass');
const MassUnit = goog.require('proto.ord.Mass.MassUnit');
const Moles = goog.require('proto.ord.Moles');
const MolesUnit = goog.require('proto.ord.Moles.MolesUnit');
const Volume = goog.require('proto.ord.Volume');
const VolumeUnit = goog.require('proto.ord.Volume.VolumeUnit');

/**
 * Initializes the selector for an Amount section.
 * @param {!jQuery} node The parent div containing the Amount section.
 */
function init(node) {
  const select = $('select.amount_units', node);
  select.append('<option value="UNSPECIFIED" selected>UNSPECIFIED</option>');
  if (select.hasClass('mass')) {
    select.append('<option disabled>-- Mass --</option>');
    for (let key in MassUnit) {
      if (key === 'UNSPECIFIED') {
        continue;
      }
      const option = $('<option>').text(key);
      option.attr('value', key);
      select.append(option);
    }
  }
  if (select.hasClass('moles')) {
    select.append('<option disabled>-- Moles --</option>');
    for (let key in MolesUnit) {
      if (key === 'UNSPECIFIED') {
        continue;
      }
      const option = $('<option>').text(key);
      option.attr('value', key);
      select.append(option);
    }
  }
  if (select.hasClass('volume')) {
    select.append('<option disabled>-- Volume --</option>');
    for (let key in VolumeUnit) {
      if (key === 'UNSPECIFIED') {
        continue;
      }
      const option = $('<option>').text(key);
      option.attr('value', key);
      select.append(option);
    }
    select.on('change', function() {
      const value = utils.getSelectorText($('.amount', node)[0]);
      if (VolumeUnit[value]) {
        $('.amount_includes_solutes_label', node).show();
        $('.amount_includes_solutes', node).show();
      } else {
        $('.amount_includes_solutes_label', node).hide();
        $('.amount_includes_solutes', node).hide();
      }
    });
  }
}

/**
 * Finds the first enum entry with the given value (see
 * https://stackoverflow.com/a/4286926).
 * @param {*} units Proto enum.
 * @param {number} value The value to find.
 * @returns {string|null}
 */
function enumFromValue(units, value) {
  for (let key in units) {
    if (units[key] === value) {
      return key;
    }
  }
  return null;
}

/**
 * Adds and populates the form's fields describing the amount of a compound.
 * @param {!jQuery} node The div corresponding to the compound whose amount
 *     fields on the form should be updated.
 * @param {?Amount} amount
 */
function load(node, amount) {
  if (!amount) {
    return;
  }
  const amountNode = $('.amount', node).first();
  const select = $('select.amount_units', amountNode);
  $('.amount_includes_solutes', node).hide();
  if (amount.hasMass()) {
    const mass = amount.getMass();
    if (mass.hasValue()) {
      $('.amount_value', node).text(mass.getValue());
    }
    if (mass.hasPrecision()) {
      $('.amount_precision', node).text(mass.getPrecision());
    }
    select.val(enumFromValue(MassUnit, mass.getUnits()));
  } else if (amount.hasMoles()) {
    const moles = amount.getMoles();
    if (moles.hasValue()) {
      $('.amount_value', node).text(moles.getValue());
    }
    if (moles.hasPrecision()) {
      $('.amount_precision', node).text(moles.getPrecision());
    }
    select.val(enumFromValue(MolesUnit, moles.getUnits()));
  } else if (amount.hasVolume()) {
    const volume = amount.getVolume();
    if (volume.hasValue()) {
      $('.amount_value', node).text(volume.getValue());
    }
    if (volume.hasPrecision()) {
      $('.amount_precision', node).text(volume.getPrecision());
    }
    select.val(enumFromValue(VolumeUnit, volume.getUnits()));
    $('.amount_includes_solutes_label', node)
        .show()
        .css('display', 'inline-block');
    $('.amount_includes_solutes', node).show().css('display', 'inline-block');
    const solutes = amount.hasVolumeIncludesSolutes() ?
        amount.getVolumeIncludesSolutes() :
        null;
    utils.setOptionalBool(
        $('.amount_includes_solutes.optional_bool', node), solutes);
  }
}

/**
 * Creates an Amount message according to the form.
 * @param {!jQuery} node The parent node for the amount fields.
 * @return {!Amount}
 */
function unload(node) {
  const amount = new Amount();
  // NOTE(kearnes): Take the closest amount section; there may be others
  // nested deeper (e.g. in ProductMeasurement fields under a ReactionProduct).
  node = $('.amount', node).first();
  const value = parseFloat($('.amount_value', node).text());
  const precision = parseFloat($('.amount_precision', node).text());
  const units = $('.amount_units', node).val();
  if (MassUnit[units]) {
    const message = new Mass();
    if (!isNaN(value)) {
      message.setValue(value);
    }
    if (!isNaN(precision)) {
      message.setPrecision(precision);
    }
    message.setUnits(MassUnit[units]);
    if (!utils.isEmptyMessage(message)) {
      amount.setMass(message);
    }
  } else if (MolesUnit[units]) {
    const message = new Moles();
    if (!isNaN(value)) {
      message.setValue(value);
    }
    if (!isNaN(precision)) {
      message.setPrecision(precision);
    }
    message.setUnits(MolesUnit[units]);
    if (!utils.isEmptyMessage(message)) {
      amount.setMoles(message);
    }
  } else if (VolumeUnit[units]) {
    const message = new Volume();
    if (!isNaN(value)) {
      message.setValue(value);
    }
    if (!isNaN(precision)) {
      message.setPrecision(precision);
    }
    message.setUnits(VolumeUnit[units]);
    if (!utils.isEmptyMessage(message)) {
      amount.setVolume(message);
    }
    const solutes = utils.getOptionalBool(
        $('.amount_includes_solutes.optional_bool', node));
    if (solutes !== null) {
      amount.setVolumeIncludesSolutes(solutes);
    }
  }
  return amount;
}
