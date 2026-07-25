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

import { render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { MemoryRouter, Route, Routes } from 'react-router-dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import ReactionCard from './ReactionCard';
import type { ReactionData, SearchResult } from '../types/search';

const reaction = (data: unknown = {}): SearchResult => ({
  reaction_id: 'ord-1',
  dataset_id: 'ord_dataset-1',
  proto: '',
  data: data as ReactionData,
});

const renderCard = (result: SearchResult, props = {}) =>
  render(
    <MemoryRouter initialEntries={['/search']}>
      <Routes>
        <Route
          path="/search"
          element={
            <ReactionCard
              reaction={result}
              {...props}
            />
          }
        />
        <Route
          path="/id/:id"
          element={<div>reaction detail page</div>}
        />
      </Routes>
    </MemoryRouter>,
  );

beforeEach(() => {
  vi.stubGlobal(
    'fetch',
    vi.fn(async () => ({ ok: true, status: 200, text: async () => '<b>A + B</b>' })),
  );
});

afterEach(() => {
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});

describe('ReactionCard', () => {
  describe('yield and conversion', () => {
    it('shows the yield measurement', async () => {
      renderCard(
        reaction({
          outcomesList: [
            {
              productsList: [
                { measurementsList: [{ type: 3, percentage: { value: 82.35 } }] },
              ],
            },
          ],
        }),
      );
      expect(await screen.findByText('Yield: 82.4%')).toBeInTheDocument();
    });

    // Only measurements of type YIELD count; the rest are purity, area, etc.
    it('ignores measurements that are not yields', async () => {
      renderCard(
        reaction({
          outcomesList: [
            {
              productsList: [
                { measurementsList: [{ type: 5, percentage: { value: 99 } }] },
              ],
            },
          ],
        }),
      );
      expect(await screen.findByText('Yield: Not listed')).toBeInTheDocument();
    });

    it('shows the conversion', async () => {
      renderCard(
        reaction({ outcomesList: [{ conversion: { value: 50, precision: 5 } }] }),
      );
      expect(await screen.findByText('Conversion: 50% ± 5')).toBeInTheDocument();
    });

    it('says so when neither was recorded', async () => {
      renderCard(reaction());
      expect(await screen.findByText('Yield: Not listed')).toBeInTheDocument();
      expect(screen.getByText('Conversion: Not listed')).toBeInTheDocument();
    });
  });

  describe('conditions', () => {
    it('joins temperature, pressure and duration', async () => {
      renderCard(
        reaction({
          conditions: {
            temperature: { setpoint: { value: 25, units: 1 } },
            pressure: { setpoint: { value: 1, units: 2 } },
          },
          outcomesList: [{ reactionTime: { value: 12, units: 1 } }],
        }),
      );
      expect(
        await screen.findByText(
          'Conditions: at 25 celsius; under 1 atmosphere; for 12 hour',
        ),
      ).toBeInTheDocument();
    });

    it('falls back to default units when the enum is unrecognized', async () => {
      renderCard(
        reaction({
          conditions: {
            temperature: { setpoint: { value: 25, units: 99 } },
            pressure: { setpoint: { value: 1, units: 99 } },
          },
          outcomesList: [{ reactionTime: { value: 12, units: 99 } }],
        }),
      );
      expect(
        await screen.findByText('Conditions: at 25°C; under 1 atm; for 12s'),
      ).toBeInTheDocument();
    });

    it('says so when none were recorded', async () => {
      renderCard(reaction());
      expect(await screen.findByText('Conditions: Not Listed')).toBeInTheDocument();
    });
  });

  describe('product identity', () => {
    it('labels the identifier with its type', async () => {
      renderCard(
        reaction({
          outcomesList: [
            { productsList: [{ identifiersList: [{ type: 2, value: 'CCO' }] }] },
          ],
        }),
      );
      expect(await screen.findByText('Product SMILES: CCO')).toBeInTheDocument();
    });

    it('omits the product row when there is no identifier', async () => {
      renderCard(reaction());
      await screen.findByText('Yield: Not listed');
      expect(screen.queryByText(/^Product /)).not.toBeInTheDocument();
    });
  });

  describe('provenance', () => {
    it('shows the uploader, date, DOI and publication link', async () => {
      renderCard(
        reaction({
          provenance: {
            doi: '10.1021/jacs.8b01523',
            publicationUrl: 'https://example.com/paper',
            recordCreated: {
              person: { name: 'Ada Lovelace', organization: 'ORD' },
              time: { value: '2021-05-04T00:00:00Z' },
            },
          },
        }),
      );

      expect(
        await screen.findByText('Uploaded by Ada Lovelace, ORD'),
      ).toBeInTheDocument();
      expect(screen.getByText('DOI: 10.1021/jacs.8b01523')).toBeInTheDocument();
      expect(screen.getByRole('link', { name: 'Publication URL' })).toHaveAttribute(
        'href',
        'https://example.com/paper',
      );
    });

    it('falls back to "Unknown" for a missing uploader', async () => {
      renderCard(reaction());
      expect(
        await screen.findByText('Uploaded by Unknown, Unknown'),
      ).toBeInTheDocument();
      expect(screen.getByText('Uploaded on Unknown')).toBeInTheDocument();
      expect(screen.getByText('DOI: Not available')).toBeInTheDocument();
    });

    it('omits the publication link when there is no URL', async () => {
      renderCard(reaction());
      await screen.findByText('DOI: Not available');
      expect(
        screen.queryByRole('link', { name: 'Publication URL' }),
      ).not.toBeInTheDocument();
    });

    it('badges a mined record', async () => {
      renderCard(reaction({ provenance: { isMined: true } }));
      expect(await screen.findByText('Mined')).toBeInTheDocument();
    });

    it('leaves an unmined record unbadged', async () => {
      renderCard(reaction({ provenance: { isMined: false } }));
      await screen.findByText('DOI: Not available');
      expect(screen.queryByText('Mined')).not.toBeInTheDocument();
    });
  });

  it('links back to the parent dataset', async () => {
    renderCard(reaction());
    expect(await screen.findByRole('link', { name: 'ord_dataset-1' })).toHaveAttribute(
      'href',
      '/search?dataset_id=ord_dataset-1&limit=100',
    );
  });

  describe('the reaction summary', () => {
    it('renders the HTML the API returns', async () => {
      const { container } = renderCard(reaction());
      await waitFor(() =>
        expect(container.querySelector('.reaction-table')?.innerHTML).toContain(
          '<b>A + B</b>',
        ),
      );
      expect(fetch).toHaveBeenCalledWith('/api/reaction_summary?reaction_id=ord-1');
    });

    // The 4xx/5xx body is an HTML error page, which would otherwise be injected
    // into every card on the results list.
    it('does not render an error page', async () => {
      const consoleError = vi.spyOn(console, 'error').mockImplementation(() => {});
      const text = vi.fn(async () => '<h1>500 Internal Server Error</h1>');
      vi.stubGlobal(
        'fetch',
        vi.fn(async () => ({ ok: false, status: 500, text })),
      );

      const { container } = renderCard(reaction());

      await waitFor(() => expect(consoleError).toHaveBeenCalled());
      expect(text).not.toHaveBeenCalled();
      expect(container.querySelector('.reaction-table')?.innerHTML).not.toContain(
        '500 Internal Server Error',
      );
    });

    it('survives a network error', async () => {
      const consoleError = vi.spyOn(console, 'error').mockImplementation(() => {});
      vi.stubGlobal(
        'fetch',
        vi.fn(async () => Promise.reject(new Error('offline'))),
      );

      renderCard(reaction());

      await waitFor(() => expect(consoleError).toHaveBeenCalled());
      expect(screen.getByText('Yield: Not listed')).toBeInTheDocument();
    });
  });

  describe('selection', () => {
    it('reports being checked and unchecked', async () => {
      const user = userEvent.setup();
      const onSelectionChange = vi.fn();
      renderCard(reaction(), { onSelectionChange, isSelected: false });

      await user.click(await screen.findByLabelText('Select reaction'));

      expect(onSelectionChange).toHaveBeenCalledWith('ord-1', true);
    });

    it('reflects the selected state it was given', async () => {
      renderCard(reaction(), { isSelected: true, onSelectionChange: vi.fn() });
      expect(await screen.findByLabelText('Select reaction')).toBeChecked();
    });

    it('hides the checkbox when the card is not selectable', async () => {
      renderCard(reaction(), { isSelectable: false });
      await screen.findByText('Yield: Not listed');
      expect(screen.queryByLabelText('Select reaction')).not.toBeInTheDocument();
    });
  });

  it('navigates to the detail page', async () => {
    const user = userEvent.setup();
    renderCard(reaction());

    await user.click(await screen.findByRole('button', { name: 'View Full Details' }));

    expect(screen.getByText('reaction detail page')).toBeInTheDocument();
  });
});
