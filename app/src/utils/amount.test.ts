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

import type { Amount } from 'ord-schema/proto/reaction_pb';
import { describe, expect, it } from 'vitest';
import { amountObj, amountStr } from './amount';

// Only the oneof field under test is populated, and some cases carry a unit
// value outside the generated enum, so the fields are cast rather than typed.
const amount = (fields: unknown): Amount.AsObject => fields as Amount.AsObject;

describe('amountObj', () => {
  it('reports no category without an amount', () => {
    expect(amountObj(undefined)).toEqual({ unitCategory: '' });
  });

  it('reports no category when no oneof field is set', () => {
    expect(amountObj(amount({}))).toEqual({ unitCategory: '' });
  });

  it('normalizes moles', () => {
    expect(
      amountObj(amount({ moles: { value: 1.5, units: 2, precision: 0 } })),
    ).toEqual({ unitAmount: 1.5, unitType: 'MILLIMOLE', unitCategory: 'moles' });
  });

  it('normalizes volume', () => {
    expect(
      amountObj(amount({ volume: { value: 10, units: 2, precision: 0 } })),
    ).toEqual({ unitAmount: 10, unitType: 'MILLILITER', unitCategory: 'volume' });
  });

  it('normalizes mass', () => {
    expect(amountObj(amount({ mass: { value: 250, units: 3, precision: 0 } }))).toEqual(
      {
        unitAmount: 250,
        unitType: 'MILLIGRAM',
        unitCategory: 'mass',
      },
    );
  });

  it('normalizes an unmeasured amount, which carries no value', () => {
    expect(amountObj(amount({ unmeasured: { type: 1, details: '' } }))).toEqual({
      unitCategory: 'unmeasured',
    });
  });

  it('leaves the unit type undefined when the units are unrecognized', () => {
    expect(amountObj(amount({ mass: { value: 1, units: 99, precision: 0 } }))).toEqual({
      unitAmount: 1,
      unitType: undefined,
      unitCategory: 'mass',
    });
  });

  it('prefers moles when more than one oneof field survived serialization', () => {
    expect(
      amountObj(
        amount({
          moles: { value: 1, units: 1, precision: 0 },
          volume: { value: 2, units: 1, precision: 0 },
        }),
      ).unitCategory,
    ).toBe('moles');
  });
});

describe('amountStr', () => {
  it('renders the value with a lowercase unit', () => {
    expect(
      amountStr({ unitAmount: 1.5, unitType: 'MILLIMOLE', unitCategory: 'moles' }),
    ).toBe('1.5 millimole');
  });

  it('rounds to three decimals', () => {
    expect(
      amountStr({ unitAmount: 1.23456, unitType: 'GRAM', unitCategory: 'mass' }),
    ).toBe('1.235 gram');
    expect(
      amountStr({ unitAmount: 0.0001, unitType: 'GRAM', unitCategory: 'mass' }),
    ).toBe('0 gram');
  });

  // A zero amount is real data; only a *missing* amount renders as empty.
  it('renders a zero amount', () => {
    expect(amountStr({ unitAmount: 0, unitType: 'GRAM', unitCategory: 'mass' })).toBe(
      '0 gram',
    );
  });

  it('renders nothing without a value or a unit', () => {
    expect(amountStr({ unitCategory: 'unmeasured' })).toBe('');
    expect(amountStr({ unitAmount: 5, unitCategory: 'mass' })).toBe('');
    expect(amountStr({ unitType: 'GRAM', unitCategory: 'mass' })).toBe('');
  });

  it('round-trips an Amount through both helpers', () => {
    expect(
      amountStr(amountObj(amount({ volume: { value: 2.5, units: 3, precision: 0 } }))),
    ).toBe('2.5 microliter');
  });
});
