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
import type { Time } from 'ord-schema/proto/reaction_pb';
import { enumName } from './enum';

/**
 * Format a Time.AsObject as "<value> <unit>(s)" matching the Vue port's
 * `outcomesUtil.formattedTime`. Returns null when there's nothing to show.
 */
export const formattedTime = (time: Time.AsObject | undefined): string | null => {
  if (!time) return null;
  const type = enumName(reaction_pb.Time.TimeUnit, time.units);
  if (!type) return null;
  // UNSPECIFIED has enum value 0; the Vue util only pluralizes the others.
  const pluralized = time.units !== 0 ? '(s)' : '';
  return `${time.value} ${type.toLowerCase()}${pluralized}`;
};
