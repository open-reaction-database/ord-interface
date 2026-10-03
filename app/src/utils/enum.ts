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

import type { DescEnum } from '@bufbuild/protobuf';

/**
 * Look up the protobuf name of an enum value, e.g. `HOUR` for a `Time_TimeUnit`.
 *
 * Takes the enum's descriptor (`Time_TimeUnitSchema`) so the name comes from the
 * schema. Proto3 enums are open, so a decoded record can carry a number the schema
 * does not declare; that returns undefined.
 */
export function enumName(
  schema: DescEnum,
  value: number | undefined,
): string | undefined {
  if (value === undefined) return undefined;
  return schema.value[value]?.name;
}
