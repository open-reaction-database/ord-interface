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

import reaction_pb from 'ord-schema';
import type {
  ElectrochemistryConditions,
  FlowConditions,
  IlluminationConditions,
  Length,
  PressureConditions,
  StirringConditions,
  Temperature,
  Pressure,
} from 'ord-schema/proto/reaction_pb';
import { enumName } from './enum';

export const tempType = (type: number | undefined): string =>
  enumName(reaction_pb.TemperatureConditions.TemperatureControl.TemperatureControlType, type) ?? '';

export const tempSetPoint = (setpoint: Temperature.AsObject | undefined): string => {
  if (!setpoint) return 'None';
  const unit = enumName(reaction_pb.Temperature.TemperatureUnit, setpoint.units);
  const precision = setpoint.precision ? ` (± ${setpoint.precision})` : '';
  return `${setpoint.value}${precision} °${unit ? unit.charAt(0) : ''}`;
};

export const pressureType = (type: number | undefined): string =>
  enumName(reaction_pb.PressureConditions.PressureControl.PressureControlType, type) ?? '';

export const pressureSetPoint = (setpoint: Pressure.AsObject | undefined): string => {
  if (!setpoint) return 'None';
  const unit = enumName(reaction_pb.Pressure.PressureUnit, setpoint.units);
  const precision = setpoint.precision ? ` (± ${setpoint.precision})` : '';
  return `${setpoint.value}${precision} ${unit ? unit.toLowerCase() : ''}`;
};

export const pressureAtmo = (atmo: PressureConditions.Atmosphere.AsObject | undefined): string => {
  const type = enumName(reaction_pb.PressureConditions.Atmosphere.AtmosphereType, atmo?.type);
  return `${type ?? ''}${atmo?.details ? `, ${atmo.details}` : ''}`;
};

export const stirType = (type: number | undefined): string =>
  enumName(reaction_pb.StirringConditions.StirringMethodType, type) ?? '';

/**
 * The Vue util mistakenly passed the whole StirringRate object and compared it
 * to numeric enum values, so "Rate" always rendered as undefined. Take the
 * type field explicitly.
 */
export const stirRate = (rate: StirringConditions.StirringRate.AsObject | undefined): string =>
  enumName(reaction_pb.StirringConditions.StirringRate.StirringRateType, rate?.type) ?? '';

export const illumType = (illum: IlluminationConditions.AsObject | undefined): string => {
  if (!illum) return '';
  const type = enumName(reaction_pb.IlluminationConditions.IlluminationType, illum.type);
  return `${type ?? ''}${illum.details ? `: ${illum.details}` : ''}`;
};

/**
 * Format a Length.AsObject like "5 mm". Returns `undefined` when there's
 * nothing to show. The Vue ConditionsView used to render the whole Length
 * object, which serialized as "[object Object]".
 */
export const lengthStr = (length: Length.AsObject | undefined): string | undefined => {
  if (!length) return undefined;
  const unit = enumName(reaction_pb.Length.LengthUnit, length.units);
  return `${length.value}${unit ? ` ${unit.toLowerCase()}` : ''}`;
};

export const electrochemType = (type: ElectrochemistryConditions.AsObject['type'] | undefined): string =>
  enumName(reaction_pb.ElectrochemistryConditions.ElectrochemistryType, type) ?? '';

export const flowType = (type: FlowConditions.AsObject['type'] | undefined): string =>
  enumName(reaction_pb.FlowConditions.FlowType, type) ?? '';
