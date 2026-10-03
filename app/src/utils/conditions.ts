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

import {
  ElectrochemistryConditions_ElectrochemistryTypeSchema,
  FlowConditions_FlowTypeSchema,
  IlluminationConditions_IlluminationTypeSchema,
  Length_LengthUnitSchema,
  PressureConditions_Atmosphere_AtmosphereTypeSchema,
  PressureConditions_PressureControl_PressureControlTypeSchema,
  Pressure_PressureUnitSchema,
  StirringConditions_StirringMethodTypeSchema,
  StirringConditions_StirringRate_StirringRateTypeSchema,
  TemperatureConditions_TemperatureControl_TemperatureControlTypeSchema,
  Temperature_TemperatureUnitSchema,
  Wavelength_WavelengthUnitSchema,
  type ElectrochemistryConditions,
  type FlowConditions,
  type IlluminationConditions,
  type Length,
  type Pressure,
  type PressureConditions_Atmosphere,
  type StirringConditions_StirringRate,
  type Temperature,
  type Wavelength,
} from '@buf/open-reaction-database_ord-schema.bufbuild_es/ord-schema/proto/reaction_pb';
import { enumName } from './enum';

// Setpoints, lengths, and wavelengths read an unset value as 0, the proto3 default.

export const tempType = (type: number | undefined): string =>
  enumName(
    TemperatureConditions_TemperatureControl_TemperatureControlTypeSchema,
    type,
  ) ?? '';

export const tempSetPoint = (setpoint: Temperature | undefined): string => {
  if (!setpoint) return 'None';
  const unit = enumName(Temperature_TemperatureUnitSchema, setpoint.units);
  const precision = setpoint.precision ? ` (± ${setpoint.precision})` : '';
  return `${setpoint.value ?? 0}${precision} °${unit ? unit.charAt(0) : ''}`;
};

export const pressureType = (type: number | undefined): string =>
  enumName(PressureConditions_PressureControl_PressureControlTypeSchema, type) ?? '';

export const pressureSetPoint = (setpoint: Pressure | undefined): string => {
  if (!setpoint) return 'None';
  const unit = enumName(Pressure_PressureUnitSchema, setpoint.units);
  const precision = setpoint.precision ? ` (± ${setpoint.precision})` : '';
  return `${setpoint.value ?? 0}${precision} ${unit ? unit.toLowerCase() : ''}`;
};

export const pressureAtmo = (
  atmo: PressureConditions_Atmosphere | undefined,
): string => {
  const type = enumName(PressureConditions_Atmosphere_AtmosphereTypeSchema, atmo?.type);
  return `${type ?? ''}${atmo?.details ? `, ${atmo.details}` : ''}`;
};

export const stirType = (type: number | undefined): string =>
  enumName(StirringConditions_StirringMethodTypeSchema, type) ?? '';

/**
 * The Vue util mistakenly passed the whole StirringRate object and compared it
 * to numeric enum values, so "Rate" always rendered as undefined. Take the
 * type field explicitly.
 */
export const stirRate = (rate: StirringConditions_StirringRate | undefined): string =>
  enumName(StirringConditions_StirringRate_StirringRateTypeSchema, rate?.type) ?? '';

export const illumType = (illum: IlluminationConditions | undefined): string => {
  if (!illum) return '';
  const type = enumName(IlluminationConditions_IlluminationTypeSchema, illum.type);
  return `${type ?? ''}${illum.details ? `: ${illum.details}` : ''}`;
};

/**
 * Format a Length like "5 millimeter". Returns `undefined` when there's
 * nothing to show. The Vue ConditionsView used to render the whole Length
 * object, which serialized as "[object Object]".
 */
export const lengthStr = (length: Length | undefined): string | undefined => {
  if (!length) return undefined;
  const unit = enumName(Length_LengthUnitSchema, length.units);
  return `${length.value ?? 0}${unit ? ` ${unit.toLowerCase()}` : ''}`;
};

/**
 * Format a Wavelength like "450 nanometer". Returns `undefined` when there's
 * nothing to show.
 */
export const wavelengthStr = (
  wavelength: Wavelength | undefined,
): string | undefined => {
  if (!wavelength) return undefined;
  const unit = enumName(Wavelength_WavelengthUnitSchema, wavelength.units);
  return `${wavelength.value ?? 0}${unit ? ` ${unit.toLowerCase()}` : ''}`;
};

export const electrochemType = (
  type: ElectrochemistryConditions['type'] | undefined,
): string =>
  enumName(ElectrochemistryConditions_ElectrochemistryTypeSchema, type) ?? '';

export const flowType = (type: FlowConditions['type'] | undefined): string =>
  enumName(FlowConditions_FlowTypeSchema, type) ?? '';
