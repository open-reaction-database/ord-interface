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

import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import OutcomesView from './OutcomesView';
import type { ReactionOutcomeData } from '../../types/search';

const YIELD = 3;
const CUSTOM = 1;
const AMOUNT = 9;
const LC_ANALYSIS = 2;

const renderOutcome = (outcome: unknown) =>
  render(<OutcomesView outcome={outcome as ReactionOutcomeData} />);

const detailFields = (container: HTMLElement): Record<string, string> => {
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

// The measurements grid is five columns wide: type, spacer, value, analysis, raw.
const measurementRows = (container: HTMLElement): string[][] => {
  const cells = [...(container.querySelector('.measurements')?.children ?? [])].map(
    cell => cell.textContent ?? '',
  );
  const rows: string[][] = [];
  for (let index = 5; index < cells.length; index += 5) {
    rows.push(cells.slice(index, index + 5));
  }
  return rows;
};

const product = (measurementsList: unknown[]) => ({
  identifiersList: [],
  measurementsList,
});

beforeEach(() => {
  // CompoundView, rendered for the selected product, asks for a drawing.
  vi.stubGlobal(
    'fetch',
    vi.fn(async () => ({ ok: true, status: 200, json: async () => '<svg></svg>' })),
  );
});

afterEach(() => {
  vi.unstubAllGlobals();
});

describe('OutcomesView', () => {
  it('renders nothing without an outcome', () => {
    const { container } = renderOutcome(undefined);
    expect(container).toBeEmptyDOMElement();
  });

  describe('details', () => {
    it('shows the reaction time and conversion', () => {
      const { container } = renderOutcome({
        reactionTime: { value: 2, units: 1, precision: 0 },
        conversion: { value: 95, precision: 2 },
      });

      expect(detailFields(container)).toEqual({
        'Reaction Time': '2 hour(s)',
        Conversion: '95% ± 2',
      });
    });

    it('omits the details section when neither was recorded', () => {
      renderOutcome({ productsList: [] });
      expect(screen.queryByText('Details')).not.toBeInTheDocument();
    });
  });

  describe('products', () => {
    it('renders a tab per product and opens the first', () => {
      const { container } = renderOutcome({
        productsList: [product([]), product([])],
      });

      const tabs = [...container.querySelectorAll('.tab')];
      expect(tabs.map(tab => tab.textContent)).toEqual(['Product 1', 'Product 2']);
      expect(tabs[0]).toHaveClass('selected');
    });

    it('switches products on click', async () => {
      const user = userEvent.setup();
      const { container } = renderOutcome({
        productsList: [
          product([{ type: YIELD, percentage: { value: 10, precision: 0 } }]),
          product([{ type: YIELD, percentage: { value: 20, precision: 0 } }]),
        ],
      });

      await user.click([...container.querySelectorAll('.tab')][1]);

      expect(measurementRows(container)[0][2]).toBe('20%');
    });
  });

  describe('measurements', () => {
    it('renders a percentage measurement', () => {
      const { container } = renderOutcome({
        productsList: [
          product([
            {
              type: YIELD,
              percentage: { value: 82.35, precision: 0 },
              analysisKey: 'lcms',
            },
          ]),
        ],
      });

      expect(measurementRows(container)).toEqual([
        ['YIELD', '', '82.4%', 'lcms', '<>'],
      ]);
    });

    it('renders an amount measurement', () => {
      const { container } = renderOutcome({
        productsList: [
          product([
            {
              type: AMOUNT,
              amount: { mass: { value: 250, units: 3, precision: 0 } },
              analysisKey: '',
            },
          ]),
        ],
      });

      expect(measurementRows(container)[0][2]).toBe('250 milligram');
    });

    it('renders a float measurement', () => {
      const { container } = renderOutcome({
        productsList: [
          product([{ type: AMOUNT, floatValue: { value: 1.5 }, analysisKey: '' }]),
        ],
      });

      expect(measurementRows(container)[0][2]).toBe('1.5');
    });

    it('renders a string measurement', () => {
      const { container } = renderOutcome({
        productsList: [
          product([{ type: AMOUNT, stringValue: 'trace', analysisKey: '' }]),
        ],
      });

      expect(measurementRows(container)[0][2]).toBe('trace');
    });

    it('leaves the value blank when the oneof is empty', () => {
      const { container } = renderOutcome({
        productsList: [product([{ type: AMOUNT, analysisKey: '' }])],
      });

      expect(measurementRows(container)[0][2]).toBe('');
    });

    it('opens the raw measurement with its type named', async () => {
      const user = userEvent.setup();
      const { container } = renderOutcome({
        productsList: [
          product([
            {
              type: YIELD,
              percentage: { value: 50, precision: 0 },
              analysisKey: 'nmr',
            },
          ]),
        ],
      });

      await user.click(container.querySelector('.measurements .button')!);

      expect(screen.getByText('Raw Data')).toBeInTheDocument();
      expect(screen.getByText(/"type": "YIELD"/)).toBeInTheDocument();
    });

    describe('custom measurements', () => {
      it('offers the recorded details', async () => {
        const user = userEvent.setup();
        renderOutcome({
          productsList: [
            product([{ type: CUSTOM, details: 'in-house assay', analysisKey: '' }]),
          ],
        });

        await user.click(screen.getByText('CUSTOM'));

        expect(screen.getByText('Custom Measurement Details')).toBeInTheDocument();
        expect(screen.getByText('in-house assay')).toBeInTheDocument();
      });

      it('points at the author when no details were recorded', async () => {
        const user = userEvent.setup();
        renderOutcome({
          productsList: [product([{ type: CUSTOM, details: '', analysisKey: '' }])],
        });

        await user.click(screen.getByText('CUSTOM'));

        expect(
          screen.getByText(/contact the author for details on this custom measurement/),
        ).toBeInTheDocument();
      });

      it('closes the details again', async () => {
        const user = userEvent.setup();
        renderOutcome({
          productsList: [
            product([{ type: CUSTOM, details: 'in-house assay', analysisKey: '' }]),
          ],
        });

        await user.click(screen.getByText('CUSTOM'));
        await user.click(screen.getByText('✕'));

        expect(
          screen.queryByText('Custom Measurement Details'),
        ).not.toBeInTheDocument();
      });
    });
  });

  describe('analyses', () => {
    it('omits the section when there are none', () => {
      renderOutcome({ productsList: [product([])], analysesMap: [] });
      expect(screen.queryByText('Analyses')).not.toBeInTheDocument();
    });

    it('names each analysis by its map key and shows the selected one', () => {
      const { container } = renderOutcome({
        productsList: [product([])],
        analysesMap: [
          ['lcms', { type: LC_ANALYSIS, details: 'method A' }],
          ['nmr', { type: 5, details: 'method B' }],
        ],
      });

      const tabs = [...container.querySelectorAll('.tab')];
      expect(tabs.map(tab => tab.textContent)).toEqual(['Product 1', 'lcms', 'nmr']);
      expect(screen.getByText('LC')).toBeInTheDocument();
      expect(screen.getByText('method A')).toBeInTheDocument();
    });

    it('switches analyses on click', async () => {
      const user = userEvent.setup();
      const { container } = renderOutcome({
        productsList: [product([])],
        analysesMap: [
          ['lcms', { type: LC_ANALYSIS, details: 'method A' }],
          ['nmr', { type: 5, details: 'method B' }],
        ],
      });

      await user.click([...container.querySelectorAll('.tab')][2]);

      expect(screen.getByText('NMR_1H')).toBeInTheDocument();
      expect(screen.getByText('method B')).toBeInTheDocument();
    });

    it('opens the raw analysis with its type named', async () => {
      const user = userEvent.setup();
      const { container } = renderOutcome({
        productsList: [product([])],
        analysesMap: [['lcms', { type: LC_ANALYSIS, details: 'method A' }]],
      });

      await user.click(container.querySelector('.details .button')!);

      expect(screen.getByText(/"type": "LC"/)).toBeInTheDocument();
    });
  });
});
