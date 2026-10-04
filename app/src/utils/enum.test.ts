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
  Mass_MassUnitSchema,
  Pressure_PressureUnitSchema,
  Time_TimeUnit,
  Time_TimeUnitSchema,
} from '@buf/open-reaction-database_ord-schema.bufbuild_es/ord-schema/proto/reaction_pb';
import { describe, expect, it } from 'vitest';
import { enumName } from './enum';

describe('enumName', () => {
  it('names a value from the enum schema', () => {
    expect(enumName(Time_TimeUnitSchema, 1)).toBe('HOUR');
    expect(enumName(Mass_MassUnitSchema, 2)).toBe('GRAM');
    expect(enumName(Pressure_PressureUnitSchema, 8)).toBe('MM_HG');
  });

  // DAY is declared second but numbered 4, so a lookup by position would miss it.
  it('looks values up by number, not by declaration order', () => {
    expect(enumName(Time_TimeUnitSchema, Time_TimeUnit.DAY)).toBe('DAY');
  });

  // Zero is a real enum value in proto3, not a missing one.
  it('resolves the zero value', () => {
    expect(enumName(Time_TimeUnitSchema, 0)).toBe('UNSPECIFIED');
  });

  it('returns undefined for a value the enum does not declare', () => {
    expect(enumName(Time_TimeUnitSchema, 99)).toBeUndefined();
    expect(enumName(Time_TimeUnitSchema, -1)).toBeUndefined();
  });

  it('returns undefined without a value', () => {
    expect(enumName(Time_TimeUnitSchema, undefined)).toBeUndefined();
  });
});
