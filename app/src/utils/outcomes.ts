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

/**
 * Format a Time message as "<value> <unit>(s)" matching the Vue port's
 * `outcomesUtil.formattedTime`. Returns null when there's nothing to show.
 */
export const formattedTime = (time: ord.ITime | null | undefined): string | null => {
  if (!time) return null;
  const type = enumName(ord.Time.TimeUnit, time.units);
  if (!type) return null;
  // UNSPECIFIED has enum value 0; the Vue util only pluralizes the others.
  const pluralized = time.units !== 0 ? '(s)' : '';
  return `${time.value ?? 0} ${type.toLowerCase()}${pluralized}`;
};

/**
 * Format a Percentage message as "X%" or "X% ± Y", rounded to one decimal.
 * Used by both ReactionCard yield/conversion and OutcomesView so the two
 * call sites stay in sync.
 */
export const formatPercentage = (
  percentage: ord.IPercentage | null | undefined,
): string => {
  if (!percentage) return '';
  const rounded = Math.round((percentage.value ?? 0) * 10) / 10;
  // Percentage.precision defaults to 0 in proto3; treat 0 as "no precision
  // recorded" rather than "± 0", matching the Vue OutcomesView's
  // `isNaN(precision)` guard. Round to one decimal: the schema stores float32,
  // so raw values carry float→double noise (4.8 arrives as 4.800000190734863).
  const precision =
    percentage.precision != null &&
    Number.isFinite(percentage.precision) &&
    percentage.precision !== 0
      ? ` ± ${Math.round(percentage.precision * 10) / 10}`
      : '';
  return `${rounded}%${precision}`;
};
