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

import { render } from '@testing-library/react';
import { describe, expect, it } from 'vitest';
import WorkupsView from './WorkupsView';
import type { ReactionWorkupData } from '../../types/search';

const SMILES = 2;
const NAME = 6;

const renderWorkup = (workup: unknown) =>
  render(<WorkupsView workup={workup as ReactionWorkupData} />);

const fields = (container: HTMLElement): Record<string, string> => {
  const details = container.querySelector('.details');
  if (!details) return {};
  const labels = [...details.querySelectorAll('.label')].map(
    el => el.textContent ?? '',
  );
  const values = [...details.querySelectorAll('.value')].map(
    el => el.textContent ?? '',
  );
  return Object.fromEntries(labels.map((label, index) => [label, values[index] ?? '']));
};

const inputs = (container: HTMLElement): string[][] => {
  const identifiers = [...container.querySelectorAll('.identifier')];
  const amounts = [...container.querySelectorAll('.amount')];
  return identifiers.map((identifier, index) => [
    identifier.textContent ?? '',
    amounts[index]?.textContent ?? '',
  ]);
};

describe('WorkupsView', () => {
  it('renders nothing without a workup', () => {
    const { container } = renderWorkup(undefined);
    expect(container).toBeEmptyDOMElement();
  });

  it('names the workup type and shows what was recorded', () => {
    const { container } = renderWorkup({
      type: 6,
      details: 'washed twice',
      duration: { value: 30, units: 2, precision: 0 },
      amount: { volume: { value: 5, units: 2, precision: 0 } },
      keepPhase: 'organic',
      targetPh: 7,
      isAutomated: true,
    });

    expect(fields(container)).toEqual({
      Type: 'EXTRACTION',
      Details: 'washed twice',
      Duration: '30 minute(s)',
      'Aliquot amount': '5 milliliter',
      'Phase kept': 'organic',
      'Target pH': '7',
      Automated: 'yes',
    });
  });

  it('shows only the type when nothing else was recorded', () => {
    const { container } = renderWorkup({ type: 7 });
    expect(fields(container)).toEqual({ Type: 'FILTRATION' });
  });

  // proto3 defaults targetPh to 0, which is indistinguishable from a real,
  // strongly acidic reading, so it is treated as unset.
  it('hides an unset target pH', () => {
    const { container } = renderWorkup({ type: 6, targetPh: 0 });
    expect(fields(container)).not.toHaveProperty('Target pH');
  });

  describe('inputs', () => {
    it('labels each component by its NAME identifier', () => {
      const { container } = renderWorkup({
        type: 6,
        input: {
          componentsList: [
            {
              identifiersList: [
                { type: SMILES, value: 'O' },
                { type: NAME, value: 'water' },
              ],
              amount: { volume: { value: 10, units: 2, precision: 0 } },
            },
          ],
        },
      });

      expect(inputs(container)).toEqual([['water', '10 milliliter']]);
    });

    // Not every compound carries a NAME, and the view must not blow up on it.
    it('falls back to the first identifier when there is no NAME', () => {
      const { container } = renderWorkup({
        type: 6,
        input: {
          componentsList: [
            {
              identifiersList: [{ type: SMILES, value: 'CCO' }],
              amount: { volume: { value: 1, units: 2, precision: 0 } },
            },
          ],
        },
      });

      expect(inputs(container)).toEqual([['CCO', '1 milliliter']]);
    });

    it('renders an empty label when the component has no identifiers', () => {
      const { container } = renderWorkup({
        type: 6,
        input: { componentsList: [{ identifiersList: [] }] },
      });

      expect(inputs(container)).toEqual([['', '']]);
    });

    it('omits the inputs section when there are no components', () => {
      const { container } = renderWorkup({ type: 6, input: { componentsList: [] } });
      expect(container.querySelector('.inputs')).toBeNull();
    });
  });
});
