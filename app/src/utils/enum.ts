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
 * Reverse-lookup the name of a protobuf enum map (e.g.
 * `ReactionSetup.ReactionEnvironment.ReactionEnvironmentType`) for a given numeric value.
 *
 * ord-schema generates these as TS interfaces (`{ UNSPECIFIED: 0; CUSTOM: 1; ... }`)
 * rather than indexable records, so this helper takes `unknown` for the map and
 * narrows internally rather than asking callers to add unsafe casts.
 */
export function enumName(
  enumMap: unknown,
  value: number | undefined,
): string | undefined {
  if (value === undefined || enumMap === null || typeof enumMap !== 'object')
    return undefined;
  for (const [key, mapped] of Object.entries(enumMap as Record<string, unknown>)) {
    if (mapped === value) return key;
  }
  return undefined;
}
