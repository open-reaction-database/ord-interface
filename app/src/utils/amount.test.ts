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

import { create, type MessageInitShape } from '@bufbuild/protobuf';
import {
  AmountSchema,
  type Mass_MassUnit,
} from '@buf/open-reaction-database_ord-schema.bufbuild_es/ord-schema/proto/reaction_pb';
import { describe, expect, it } from 'vitest';
import { amountObj, amountStr } from './amount';

const amount = (init: MessageInitShape<typeof AmountSchema>) =>
  create(AmountSchema, init);

describe('amountObj', () => {
  it('reports no category without an amount', () => {
    expect(amountObj(undefined)).toEqual({ unitCategory: '' });
  });

  it('reports no category when no oneof field is set', () => {
    expect(amountObj(amount({}))).toEqual({ unitCategory: '' });
  });

  it('normalizes moles', () => {
    expect(
      amountObj(amount({ kind: { case: 'moles', value: { value: 1.5, units: 2 } } })),
    ).toEqual({ unitAmount: 1.5, unitType: 'MILLIMOLE', unitCategory: 'moles' });
  });

  it('normalizes volume', () => {
    expect(
      amountObj(amount({ kind: { case: 'volume', value: { value: 10, units: 2 } } })),
    ).toEqual({ unitAmount: 10, unitType: 'MILLILITER', unitCategory: 'volume' });
  });

  it('normalizes mass', () => {
    expect(
      amountObj(amount({ kind: { case: 'mass', value: { value: 250, units: 3 } } })),
    ).toEqual({
      unitAmount: 250,
      unitType: 'MILLIGRAM',
      unitCategory: 'mass',
    });
  });

  it('normalizes an unmeasured amount, which carries no value', () => {
    expect(
      amountObj(amount({ kind: { case: 'unmeasured', value: { type: 1 } } })),
    ).toEqual({
      unitCategory: 'unmeasured',
    });
  });

  it('reads an unset value as zero', () => {
    expect(amountObj(amount({ kind: { case: 'mass', value: { units: 2 } } }))).toEqual({
      unitAmount: 0,
      unitType: 'GRAM',
      unitCategory: 'mass',
    });
  });

  // Proto3 enums are open, so a decoded record can carry an undeclared unit.
  it('leaves the unit type undefined when the units are unrecognized', () => {
    expect(
      amountObj(
        amount({
          kind: { case: 'mass', value: { value: 1, units: 99 as Mass_MassUnit } },
        }),
      ),
    ).toEqual({
      unitAmount: 1,
      unitType: undefined,
      unitCategory: 'mass',
    });
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
      amountStr(
        amountObj(
          amount({ kind: { case: 'volume', value: { value: 2.5, units: 3 } } }),
        ),
      ),
    ).toBe('2.5 microliter');
  });
});
