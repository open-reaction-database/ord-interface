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
  Moles_MolesUnitSchema,
  Volume_VolumeUnitSchema,
  type Amount,
} from '@buf/open-reaction-database_ord-schema.bufbuild_es/ord-schema/proto/reaction_pb';
import { enumName } from './enum';

export type AmountCategory = 'moles' | 'volume' | 'mass' | 'unmeasured' | '';

export interface AmountObj {
  unitAmount?: number;
  unitType?: string;
  unitCategory: AmountCategory;
}

/**
 * Normalize a Compound.amount oneof into a flat { unitAmount, unitType, unitCategory }
 * triple so render code doesn't have to branch on which oneof field is populated.
 *
 * An unset value reads as 0, the proto3 default.
 */
export const amountObj = (amount: Amount | undefined): AmountObj => {
  switch (amount?.kind.case) {
    case 'moles':
      return {
        unitAmount: amount.kind.value.value ?? 0,
        unitType: enumName(Moles_MolesUnitSchema, amount.kind.value.units),
        unitCategory: 'moles',
      };
    case 'volume':
      return {
        unitAmount: amount.kind.value.value ?? 0,
        unitType: enumName(Volume_VolumeUnitSchema, amount.kind.value.units),
        unitCategory: 'volume',
      };
    case 'mass':
      return {
        unitAmount: amount.kind.value.value ?? 0,
        unitType: enumName(Mass_MassUnitSchema, amount.kind.value.units),
        unitCategory: 'mass',
      };
    case 'unmeasured':
      return { unitCategory: 'unmeasured' };
    default:
      return { unitCategory: '' };
  }
};

export const amountStr = (obj: AmountObj): string => {
  if (obj.unitAmount === undefined || !obj.unitType) return '';
  return `${Math.round(obj.unitAmount * 1000) / 1000} ${obj.unitType.toLowerCase()}`;
};
