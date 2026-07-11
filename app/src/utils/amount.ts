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

export type AmountCategory = 'moles' | 'volume' | 'mass' | 'unmeasured' | '';

export interface AmountObj {
  unitAmount?: number;
  unitType?: string;
  unitCategory: AmountCategory;
}

/**
 * Normalize a Compound.amount oneof into a flat { unitAmount, unitType,
 * unitCategory } triple so render code doesn't have to branch on which oneof
 * field is populated.
 */
export const amountObj = (amount: ord.IAmount | null | undefined): AmountObj => {
  if (!amount) return { unitCategory: '' };
  if (amount.moles) {
    return {
      unitAmount: amount.moles.value ?? undefined,
      unitType: enumName(ord.Moles.MolesUnit, amount.moles.units),
      unitCategory: 'moles',
    };
  }
  if (amount.volume) {
    return {
      unitAmount: amount.volume.value ?? undefined,
      unitType: enumName(ord.Volume.VolumeUnit, amount.volume.units),
      unitCategory: 'volume',
    };
  }
  if (amount.mass) {
    return {
      unitAmount: amount.mass.value ?? undefined,
      unitType: enumName(ord.Mass.MassUnit, amount.mass.units),
      unitCategory: 'mass',
    };
  }
  if (amount.unmeasured) {
    return { unitCategory: 'unmeasured' };
  }
  return { unitCategory: '' };
};

export const amountStr = (obj: AmountObj): string => {
  if (obj.unitAmount === undefined || !obj.unitType) return '';
  return `${Math.round(obj.unitAmount * 1000) / 1000} ${obj.unitType.toLowerCase()}`;
};
