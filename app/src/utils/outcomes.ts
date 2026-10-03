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
  Time_TimeUnitSchema,
  type Percentage,
  type Time,
} from '@buf/open-reaction-database_ord-schema.bufbuild_es/ord-schema/proto/reaction_pb';
import { enumName } from './enum';

/**
 * Format a Time as "<value> <unit>(s)" matching the Vue port's
 * `outcomesUtil.formattedTime`. Returns null when there's nothing to show.
 *
 * An unset value reads as 0, the proto3 default.
 */
export const formattedTime = (time: Time | undefined): string | null => {
  if (!time) return null;
  const type = enumName(Time_TimeUnitSchema, time.units);
  if (!type) return null;
  // UNSPECIFIED has enum value 0; the Vue util only pluralizes the others.
  const pluralized = time.units !== 0 ? '(s)' : '';
  return `${time.value ?? 0} ${type.toLowerCase()}${pluralized}`;
};

/**
 * Format a Percentage as "X%" or "X% ± Y", rounded to one decimal.
 * Used by both ReactionCard yield/conversion and OutcomesView so the two
 * call sites stay in sync. An unset value reads as 0, the proto3 default.
 */
export const formatPercentage = (percentage: Percentage | undefined): string => {
  if (!percentage) return '';
  const rounded = Math.round((percentage.value ?? 0) * 10) / 10;
  // A precision of 0 renders like an unset one rather than as "± 0", matching
  // the Vue OutcomesView's `isNaN(precision)` guard.
  const precision =
    Number.isFinite(percentage.precision) && percentage.precision !== 0
      ? ` ± ${percentage.precision}`
      : '';
  return `${rounded}%${precision}`;
};
