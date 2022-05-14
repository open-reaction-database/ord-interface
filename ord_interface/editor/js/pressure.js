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

goog.module('ord.pressure');
goog.module.declareLegacyNamespace();

const asserts = goog.require('goog.asserts');

const utils = goog.require('ord.utils');

const Pressure = goog.require('proto.ord.Pressure');
const PressureConditions = goog.require('proto.ord.PressureConditions');
const Atmosphere = goog.require('proto.ord.PressureConditions.Atmosphere');
const AtmosphereType = goog.require('proto.ord.PressureConditions.Atmosphere.AtmosphereType');
const Measurement = goog.require('proto.ord.PressureConditions.Measurement');
const MeasurementType = goog.require('proto.ord.PressureConditions.Measurement.MeasurementType');
const PressureControl = goog.require('proto.ord.PressureConditions.PressureControl');
const PressureControlType = goog.require('proto.ord.PressureConditions.PressureControl.PressureControlType');
const Time = goog.require('proto.ord.Time');

exports = {
  load,
  unload,
  addMeasurement,
  validatePressure
};


/**
 * Adds and populates the pressure conditions section in the form.
 * @param {!PressureConditions} pressure
 */
function load(pressure) {
  const control = pressure.getControl();
  if (control) {
    utils.setSelector($('#pressure_control_type'), control.getType());
    $('#pressure_control_details').text(control.getDetails());
  }
  const measurements = pressure.getMeasurementsList();
  measurements.forEach(function(measurement) {
    const node = addMeasurement();
    loadMeasurement(measurement, node);
  });
  const setpoint = pressure.getSetpoint();
  utils.writeMetric('#pressure_setpoint', setpoint);

  const atmosphere = pressure.getAtmosphere();
  if (atmosphere) {
    utils.setSelector($('#pressure_atmosphere_type'), atmosphere.getType());
    $('#pressure_atmosphere_details').text(atmosphere.getDetails());
  }
}

/**
 * Adds and populates a pressure measurement section in the form.
 * @param {!PressureConditions.Measurement} measurement
 * @param {!jQuery} node The target div.
 */
function loadMeasurement(measurement, node) {
  const type = measurement.getType();
  utils.setSelector($('.pressure_measurement_type', node), type);
  $('.pressure_measurement_details', node).text(measurement.getDetails());

  const pressure = measurement.getPressure();
  utils.writeMetric('.pressure_measurement_pressure', pressure, node);

  const time = measurement.getTime();
  utils.writeMetric('.pressure_measurement_time', time, node);
}

/**
 * Fetches pressure conditions from the form.
 * @return {!PressureConditions}
 */
function unload() {
  const pressure = new PressureConditions();

  const control = new PressureControl();
  const controlType = utils.getSelectorText($('#pressure_control_type')[0]);
  control.setType(PressureControlType[controlType]);
  control.setDetails(
      asserts.assertString($('#pressure_control_details').text()));
  if (!utils.isEmptyMessage(control)) {
    pressure.setControl(control);
  }

  const setpoint = utils.readMetric('#pressure_setpoint', new Pressure());
  if (!utils.isEmptyMessage(setpoint)) {
    pressure.setSetpoint(setpoint);
  }

  const atmosphere = new Atmosphere();
  const atmosphereType =
      utils.getSelectorText($('#pressure_atmosphere_type')[0]);
  atmosphere.setType(AtmosphereType[atmosphereType]);
  atmosphere.setDetails(
      asserts.assertString($('#pressure_atmosphere_details').text()));
  if (!utils.isEmptyMessage(atmosphere)) {
    pressure.setAtmosphere(atmosphere);
  }

  const measurements = [];
  $('.pressure_measurement').each(function(index, node) {
    node = $(node);
    if (!utils.isTemplateOrUndoBuffer(node)) {
      const measurement = unloadMeasurement(node);
      if (!utils.isEmptyMessage(measurement)) {
        measurements.push(measurement);
      }
    }
  });
  pressure.setMeasurementsList(measurements);

  return pressure;
}

/**
 * Fetches a pressure measurement from the form.
 * @param {!jQuery} node The div of the measurement to fetch.
 * @return {!PressureConditions.Measurement}
 */
function unloadMeasurement(node) {
  const measurement = new Measurement();
  const measurementType =
      utils.getSelectorText($('.pressure_measurement_type', node)[0]);
  measurement.setType(MeasurementType[measurementType]);
  measurement.setDetails(
      asserts.assertString($('.pressure_measurement_details', node).text()));
  const pressure =
      utils.readMetric('.pressure_measurement_pressure', new Pressure(), node);
  if (!utils.isEmptyMessage(pressure)) {
    measurement.setPressure(pressure);
  }
  const time = utils.readMetric('.pressure_measurement_time', new Time(), node);
  if (!utils.isEmptyMessage(time)) {
    measurement.setTime(time);
  }

  return measurement;
}

/**
 * Adds a pressure measurement section to the form.
 * @return {!jQuery} The node of the new measurement div.
 */
function addMeasurement() {
  return utils.addSlowly(
      '#pressure_measurement_template', $('#pressure_measurements'));
}

/**
 * Validates pressure conditions defined in the form.
 * @param {!jQuery} node The node containing the pressure conditions div.
 * @param {?jQuery=} validateNode The target div for validation results.
 */
function validatePressure(node, validateNode = null) {
  const pressure = unload();
  utils.validate(pressure, 'PressureConditions', node, validateNode);
}
