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

import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { MemoryRouter, useLocation } from 'react-router-dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import MainSearch from './MainSearch';

let currentSearch = '';

const LocationSpy = () => {
  currentSearch = useLocation().search;
  return null;
};

const renderSearch = (search = '') => {
  const client = new QueryClient({
    defaultOptions: { queries: { retry: false, gcTime: 0 } },
  });
  currentSearch = search;
  return render(
    <QueryClientProvider client={client}>
      <MemoryRouter initialEntries={[`/search${search}`]}>
        <LocationSpy />
        <MainSearch />
      </MemoryRouter>
    </QueryClientProvider>,
  );
};

// The params of the URL the component navigated to, as sorted "key=value" pairs
// so assertions do not depend on the order they were appended in.
const navigatedParams = (): string[] =>
  [...new URLSearchParams(currentSearch).entries()]
    .map(([key, value]) => `${key}=${value}`)
    .sort();

const searchButton = () => screen.getByRole('button', { name: 'Search' });

beforeEach(() => {
  // Never resolves: the submit/poll protocol is exercised in the useSearchTask
  // tests, and here it just has to stay out of the way.
  vi.stubGlobal(
    'fetch',
    vi.fn(() => new Promise<Response>(() => {})),
  );
});

afterEach(() => {
  vi.unstubAllGlobals();
});

describe('MainSearch', () => {
  it('prompts for criteria when the URL carries none', () => {
    renderSearch();
    expect(screen.getByText(/Enter search criteria/)).toBeInTheDocument();
  });

  // `limit` alone is not a search; it would otherwise fire an empty query on
  // every visit.
  it('does not treat the result limit as search criteria', () => {
    renderSearch('?limit=100');
    expect(screen.getByText(/Enter search criteria/)).toBeInTheDocument();
    expect(fetch).not.toHaveBeenCalled();
  });

  it('runs a search when the URL carries criteria', () => {
    renderSearch('?dataset_id=ord_dataset-1');
    expect(fetch).toHaveBeenCalledWith(
      '/api/submit_query?dataset_id=ord_dataset-1',
      undefined,
    );
  });

  it('encodes components as JSON with the shared match mode', async () => {
    const user = userEvent.setup();
    renderSearch(
      `?component=${encodeURIComponent(
        JSON.stringify({ pattern: 'CCO', target: 'input', mode: 'substructure' }),
      )}`,
    );

    await user.click(searchButton());

    await waitFor(() =>
      expect(new URLSearchParams(currentSearch).getAll('component')).toEqual([
        JSON.stringify({ pattern: 'CCO', target: 'input', mode: 'substructure' }),
      ]),
    );
  });

  it('round-trips a full set of criteria through the URL', async () => {
    const user = userEvent.setup();
    const component = JSON.stringify({
      pattern: 'CCO',
      target: 'input',
      mode: 'similar',
    });
    const initial = new URLSearchParams({
      component,
      use_stereochemistry: 'true',
      similarity: '0.75',
      dataset_id: 'ord_dataset-1',
      doi: '10.1000/x',
      reaction_id: 'ord-abc',
      reaction_smarts: '[C:1]>>[C:1]O',
      min_yield: '20',
      max_yield: '80',
      min_conversion: '10',
      max_conversion: '90',
      limit: '25',
    });
    renderSearch(`?${initial.toString()}`);
    const before = navigatedParams();

    await user.click(searchButton());

    await waitFor(() => expect(navigatedParams()).toEqual(before));
  });

  // The backend applies the full range by default, so sending it back would
  // only bloat the shareable link.
  it('omits an unnarrowed yield and conversion range', async () => {
    const user = userEvent.setup();
    renderSearch('?dataset_id=ord_dataset-1');

    await user.click(searchButton());

    await waitFor(() =>
      expect(navigatedParams()).toEqual(['dataset_id=ord_dataset-1', 'limit=100']),
    );
  });

  // Stereochemistry and similarity only mean something alongside a component.
  it('omits the component settings when there are no components', async () => {
    const user = userEvent.setup();
    renderSearch('?reaction_id=ord-abc');

    await user.click(searchButton());

    await waitFor(() =>
      expect(navigatedParams()).toEqual(['limit=100', 'reaction_id=ord-abc']),
    );
  });

  it('reports a failed search', async () => {
    vi.stubGlobal(
      'fetch',
      vi.fn().mockResolvedValue({ ok: false, status: 500 } as Response),
    );
    renderSearch('?dataset_id=ord_dataset-1');

    expect(
      await screen.findByText(/Search failed: submit_query failed \(HTTP 500\)/),
    ).toBeInTheDocument();
  });

  it('reports an empty result set', async () => {
    vi.stubGlobal(
      'fetch',
      vi.fn(async (url: string) =>
        url.startsWith('/api/submit_query')
          ? ({ ok: true, status: 200, json: async () => 'task-1' } as Response)
          : ({ ok: true, status: 200, json: async () => [] } as Response),
      ),
    );
    renderSearch('?dataset_id=ord_dataset-1');

    expect(await screen.findByText(/No results/)).toBeInTheDocument();
  });
});
