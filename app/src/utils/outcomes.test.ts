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

import { create } from '@bufbuild/protobuf';
import {
  PercentageSchema,
  TimeSchema,
} from '@buf/open-reaction-database_ord-schema.bufbuild_es/ord-schema/proto/reaction_pb';
import { describe, expect, it } from 'vitest';
import { formatPercentage, formattedTime } from './outcomes';

// Proto3 enums are open, so a decoded record can carry an undeclared unit; the
// unrecognized-unit cases use 99 for one.
const time = (value: number | undefined, units: number) =>
  create(TimeSchema, { value, units });

const percentage = (value: number | undefined, precision?: number) =>
  create(PercentageSchema, { value, precision });

describe('formattedTime', () => {
  it('renders nothing without a time', () => {
    expect(formattedTime(undefined)).toBeNull();
  });

  it('renders nothing for an unrecognized unit', () => {
    expect(formattedTime(time(3, 99))).toBeNull();
  });

  it('pluralizes the unit', () => {
    expect(formattedTime(time(3, 1))).toBe('3 hour(s)');
    expect(formattedTime(time(30, 2))).toBe('30 minute(s)');
    expect(formattedTime(time(1, 4))).toBe('1 day(s)');
  });

  // "unspecified(s)" would read as nonsense, so the zero value stays singular.
  it('leaves the unspecified unit unpluralized', () => {
    expect(formattedTime(time(3, 0))).toBe('3 unspecified');
  });

  it('reads an unset value as zero', () => {
    expect(formattedTime(time(undefined, 1))).toBe('0 hour(s)');
  });
});

describe('formatPercentage', () => {
  it('renders nothing without a percentage', () => {
    expect(formatPercentage(undefined)).toBe('');
  });

  it('rounds to one decimal', () => {
    expect(formatPercentage(percentage(82.34, 0))).toBe('82.3%');
    expect(formatPercentage(percentage(82.35, 0))).toBe('82.4%');
    expect(formatPercentage(percentage(100, 0))).toBe('100%');
  });

  it('renders a zero yield rather than an empty string', () => {
    expect(formatPercentage(percentage(0, 0))).toBe('0%');
  });

  it('reads an unset value as zero', () => {
    expect(formatPercentage(percentage(undefined))).toBe('0%');
  });

  it('appends the precision when one was recorded', () => {
    expect(formatPercentage(percentage(75, 5))).toBe('75% ± 5');
  });

  // "± 0" carries no information, so a zero precision renders like an unset one.
  it('treats an unset, zero or non-finite precision as unrecorded', () => {
    expect(formatPercentage(percentage(75))).toBe('75%');
    expect(formatPercentage(percentage(75, 0))).toBe('75%');
    expect(formatPercentage(percentage(75, NaN))).toBe('75%');
    expect(formatPercentage(percentage(75, Infinity))).toBe('75%');
  });
});
