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
import { describe, expect, it } from 'vitest';
import { enumName } from './enum';

describe('enumName', () => {
  it('reverse-looks-up a generated protobuf enum map', () => {
    expect(enumName(reaction_pb.Time.TimeUnit, 1)).toBe('HOUR');
    expect(enumName(reaction_pb.Mass.MassUnit, 2)).toBe('GRAM');
    expect(enumName(reaction_pb.Pressure.PressureUnit, 8)).toBe('MM_HG');
  });

  // Zero is a real enum value in proto3, not a missing one.
  it('resolves the zero value', () => {
    expect(enumName(reaction_pb.Time.TimeUnit, 0)).toBe('UNSPECIFIED');
  });

  it('returns undefined for a value the map does not contain', () => {
    expect(enumName(reaction_pb.Time.TimeUnit, 99)).toBeUndefined();
    expect(enumName(reaction_pb.Time.TimeUnit, -1)).toBeUndefined();
  });

  it('returns undefined without a value', () => {
    expect(enumName(reaction_pb.Time.TimeUnit, undefined)).toBeUndefined();
  });

  it('returns undefined for anything that is not an enum map', () => {
    expect(enumName(null, 0)).toBeUndefined();
    expect(enumName(undefined, 0)).toBeUndefined();
    expect(enumName('HOUR', 0)).toBeUndefined();
    expect(enumName(1, 1)).toBeUndefined();
  });

  // Object.entries would happily walk an array's indices, and the numeric
  // string keys would render as garbage unit names.
  it('matches on value, not on key type', () => {
    expect(enumName({ A: 1, B: '1' }, 1)).toBe('A');
    expect(enumName({ A: '1' }, 1)).toBeUndefined();
  });
});
