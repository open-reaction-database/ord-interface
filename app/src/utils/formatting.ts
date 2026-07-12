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

/**
 * Value/precision/unit display helpers ported from ord-app
 * (renderValuePrecisionUnit + useTextFormatting + formatting.json), adapted to
 * take raw ord-schema-protobufjs messages.
 */

import { ord } from 'ord-schema-protobufjs';
import { enumName } from './enum';

/** Unit-symbol dictionary (ord-app's formatting.json, verbatim). */
const UNIT_SYMBOLS: Record<string, string> = {
  KILOGRAM: 'kg',
  GRAM: 'g',
  MILLIGRAM: 'mg',
  MICROGRAM: 'µg',
  MOLE: 'mol',
  MILLIMOLE: 'mmol',
  MICROMOLE: 'µmol',
  NANOMOLE: 'nmol',
  LITER: 'L',
  MILLILITER: 'mL',
  MICROLITER: 'µL',
  NANOLITER: 'nL',
  HOUR: 'h',
};

export function getFormattedValue(value: string): string {
  return value in UNIT_SYMBOLS ? UNIT_SYMBOLS[value] : value;
}

interface ValuePrecisionUnit {
  value?: number | null;
  precision?: number | null;
  units: string;
}

/** Formats "value ± precision unit", matching ord-app's renderValuePrecisionUnit. */
export function renderValuePrecisionUnit(
  valuePrecision: ValuePrecisionUnit | Omit<ValuePrecisionUnit, 'units'>,
): string {
  const { value, precision } = valuePrecision;
  const units = 'units' in valuePrecision ? valuePrecision.units : '';
  const precisionString = precision ? `± ${round(precision)}` : '';
  const formattedUnits = units ? getFormattedValue(units) : '';
  return [value == null ? '' : round(value), precisionString, formattedUnits]
    .filter(item => item !== '')
    .join(' ');
}

// The schema stores float32; trim float->double noise for display.
function round(value: number): number {
  return Math.round(value * 1e4) / 1e4;
}

interface MessageWithUnits {
  value?: number | null;
  precision?: number | null;
  units?: number | null;
}

/**
 * Formats a {value, precision, units} protobuf message, resolving the units
 * enum to its key (then to a symbol where one exists).
 */
export function vpu(
  message: MessageWithUnits | null | undefined,
  unitsEnum: unknown,
): string {
  if (!message) return '';
  return renderValuePrecisionUnit({
    value: message.value,
    precision: message.precision,
    units: enumName(unitsEnum, message.units) ?? '',
  });
}

/** Formats an Amount oneof with its populated branch's units. */
export function formatAmount(amount: ord.IAmount | null | undefined): string {
  if (!amount) return '';
  if (amount.moles) return vpu(amount.moles, ord.Moles.MolesUnit);
  if (amount.volume) return vpu(amount.volume, ord.Volume.VolumeUnit);
  if (amount.mass) return vpu(amount.mass, ord.Mass.MassUnit);
  if (amount.unmeasured) {
    return (
      enumName(
        ord.UnmeasuredAmount.UnmeasuredAmountType,
        amount.unmeasured.type,
      )?.replaceAll('_', ' ') ?? 'unmeasured'
    );
  }
  return '';
}

/**
 * Display value for proto3 optional bools, matching ord-app's ReactionBoolean
 * ('Unspecified' | 'True' | 'False').
 */
export function reactionBoolean(value: boolean | null | undefined): string {
  if (value == null) return 'Unspecified';
  return value ? 'True' : 'False';
}

/** Formats a DateTime message like ord-app's formatDateToDisplay (21.07.2021 03:10 pm). */
export function formatDateToDisplay(
  dateTime: ord.IDateTime | null | undefined,
): string {
  const raw = dateTime?.value;
  if (!raw) return '';
  const date = new Date(raw);
  if (Number.isNaN(date.getTime())) return raw;
  const pad = (n: number) => String(n).padStart(2, '0');
  const hours24 = date.getHours();
  const hours12 = hours24 % 12 === 0 ? 12 : hours24 % 12;
  const meridiem = hours24 < 12 ? 'am' : 'pm';
  return (
    `${pad(date.getDate())}.${pad(date.getMonth() + 1)}.${date.getFullYear()} ` +
    `${pad(hours12)}:${pad(date.getMinutes())} ${meridiem}`
  );
}
