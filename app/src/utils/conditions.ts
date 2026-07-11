/**
 * Copyright 2026 Open Reaction Database Project Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import { ord } from 'ord-schema-protobufjs';
import { enumName } from './enum';

export const tempType = (type: number | null | undefined): string =>
  enumName(ord.TemperatureConditions.TemperatureControl.TemperatureControlType, type) ??
  '';

export const tempSetPoint = (setpoint: ord.ITemperature | null | undefined): string => {
  if (!setpoint) return 'None';
  const unit = enumName(ord.Temperature.TemperatureUnit, setpoint.units);
  const precision = setpoint.precision ? ` (± ${setpoint.precision})` : '';
  return `${setpoint.value ?? 0}${precision} °${unit ? unit.charAt(0) : ''}`;
};

export const pressureType = (type: number | null | undefined): string =>
  enumName(ord.PressureConditions.PressureControl.PressureControlType, type) ?? '';

export const pressureSetPoint = (
  setpoint: ord.IPressure | null | undefined,
): string => {
  if (!setpoint) return 'None';
  const unit = enumName(ord.Pressure.PressureUnit, setpoint.units);
  const precision = setpoint.precision ? ` (± ${setpoint.precision})` : '';
  return `${setpoint.value ?? 0}${precision} ${unit ? unit.toLowerCase() : ''}`;
};

export const pressureAtmo = (
  atmo: ord.PressureConditions.IAtmosphere | null | undefined,
): string => {
  const type = enumName(ord.PressureConditions.Atmosphere.AtmosphereType, atmo?.type);
  return `${type ?? ''}${atmo?.details ? `, ${atmo.details}` : ''}`;
};

export const stirType = (type: number | null | undefined): string =>
  enumName(ord.StirringConditions.StirringMethodType, type) ?? '';

/**
 * The Vue util mistakenly passed the whole StirringRate object and compared it
 * to numeric enum values, so "Rate" always rendered as undefined. Take the
 * type field explicitly.
 */
export const stirRate = (
  rate: ord.StirringConditions.IStirringRate | null | undefined,
): string =>
  enumName(ord.StirringConditions.StirringRate.StirringRateType, rate?.type) ?? '';

export const illumType = (
  illum: ord.IIlluminationConditions | null | undefined,
): string => {
  if (!illum) return '';
  const type = enumName(ord.IlluminationConditions.IlluminationType, illum.type);
  return `${type ?? ''}${illum.details ? `: ${illum.details}` : ''}`;
};

/**
 * Format a Length message like "5 mm". Returns `undefined` when there's
 * nothing to show. The Vue ConditionsView used to render the whole Length
 * object, which serialized as "[object Object]".
 */
export const lengthStr = (
  length: ord.ILength | null | undefined,
): string | undefined => {
  if (!length) return undefined;
  const unit = enumName(ord.Length.LengthUnit, length.units);
  return `${length.value ?? 0}${unit ? ` ${unit.toLowerCase()}` : ''}`;
};

/**
 * Format a Wavelength message like "450 nanometer". Returns `undefined` when
 * there's nothing to show. Wavelength has its own unit enum, so it can't be
 * formatted with lengthStr without mislabeling the units.
 */
export const wavelengthStr = (
  wavelength: ord.IWavelength | null | undefined,
): string | undefined => {
  if (!wavelength) return undefined;
  const unit = enumName(ord.Wavelength.WavelengthUnit, wavelength.units);
  return `${wavelength.value ?? 0}${unit ? ` ${unit.toLowerCase()}` : ''}`;
};

export const electrochemType = (type: number | null | undefined): string =>
  enumName(ord.ElectrochemistryConditions.ElectrochemistryType, type) ?? '';

export const flowType = (type: number | null | undefined): string =>
  enumName(ord.FlowConditions.FlowType, type) ?? '';
