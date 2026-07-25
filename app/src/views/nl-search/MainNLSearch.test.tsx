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
import reaction_pb from 'ord-schema';
import { MemoryRouter, useLocation } from 'react-router-dom';
import { afterEach, describe, expect, it, vi } from 'vitest';
import MainNLSearch from './MainNLSearch';
import type { NLQueryResponse } from '../../types/search';

const encodedReaction = (reactionId: string): string => {
  const reaction = new reaction_pb.Reaction();
  reaction.setReactionId(reactionId);
  return btoa(String.fromCharCode(...reaction.serializeBinary()));
};

const response = (overrides: Partial<NLQueryResponse> = {}): NLQueryResponse => ({
  query: 'reactions using benzene',
  interpretation: {
    components: [{ identifier: 'benzene', target: 'INPUT', mode: 'EXACT' }],
  },
  resolved_components: [
    {
      identifier: 'benzene',
      smiles: 'c1ccccc1',
      resolver: 'pubchem',
      target: 'INPUT',
      mode: 'EXACT',
    },
  ],
  query_components: [],
  results: [{ reaction_id: 'ord-1', proto: encodedReaction('ord-1') }],
  dry_run: false,
  ...overrides,
});

let currentSearch = '';

const LocationSpy = () => {
  currentSearch = useLocation().search;
  return null;
};

const stubApi = (body: NLQueryResponse, status = 200) => {
  const fetchMock = vi.fn(async (url: string) => {
    if (url.startsWith('/api/nl_query')) {
      return { ok: status < 400, status, json: async () => body };
    }
    // Each result card asks for its own summary.
    return { ok: true, status: 200, text: async () => '' };
  });
  vi.stubGlobal('fetch', fetchMock);
  return fetchMock;
};

const renderAsk = (search = '') => {
  const client = new QueryClient({
    defaultOptions: { queries: { retry: false, gcTime: 0 } },
  });
  currentSearch = search;
  return render(
    <QueryClientProvider client={client}>
      <MemoryRouter initialEntries={[`/ask${search}`]}>
        <LocationSpy />
        <MainNLSearch />
      </MemoryRouter>
    </QueryClientProvider>,
  );
};

const queryBox = () => screen.getByRole('textbox');

afterEach(() => {
  vi.unstubAllGlobals();
});

describe('MainNLSearch', () => {
  it('flags the feature as in development', () => {
    stubApi(response());
    renderAsk();
    expect(screen.getByRole('status')).toHaveTextContent(/in development/);
  });

  it('runs nothing until a question is submitted', () => {
    const fetchMock = stubApi(response());
    renderAsk();
    expect(fetchMock).not.toHaveBeenCalled();
  });

  it('puts the submitted question in the URL so it can be shared', async () => {
    const user = userEvent.setup();
    stubApi(response());
    renderAsk();

    await user.type(queryBox(), 'reactions using benzene');
    await user.click(screen.getByRole('button', { name: 'Search' }));

    await waitFor(() =>
      expect(new URLSearchParams(currentSearch).get('q')).toBe(
        'reactions using benzene',
      ),
    );
  });

  it('trims the question before submitting', async () => {
    const user = userEvent.setup();
    stubApi(response());
    renderAsk();

    await user.type(queryBox(), '   benzene   ');
    await user.click(screen.getByRole('button', { name: 'Search' }));

    await waitFor(() =>
      expect(new URLSearchParams(currentSearch).get('q')).toBe('benzene'),
    );
  });

  it('submits an example question on click', async () => {
    const user = userEvent.setup();
    stubApi(response());
    renderAsk();

    await user.click(
      screen.getByRole('button', { name: 'reactions for synthesizing ibuprofen' }),
    );

    await waitFor(() =>
      expect(new URLSearchParams(currentSearch).get('q')).toBe(
        'reactions for synthesizing ibuprofen',
      ),
    );
  });

  it('fills the box from a question already in the URL', () => {
    stubApi(response());
    renderAsk('?q=reactions%20using%20benzene');
    expect(queryBox()).toHaveValue('reactions using benzene');
  });

  describe('the interpretation panel', () => {
    it('reports each component with its role and match mode', async () => {
      stubApi(response());
      renderAsk('?q=reactions%20using%20benzene');

      expect(await screen.findByText('reactant/reagent')).toBeInTheDocument();
      expect(screen.getByText('benzene', { selector: 'strong' })).toBeInTheDocument();
      expect(screen.getByText('(exact)')).toBeInTheDocument();
    });

    it('labels an output component as a product', async () => {
      stubApi(
        response({
          interpretation: {
            components: [
              { identifier: 'ibuprofen', target: 'OUTPUT', mode: 'SIMILAR' },
            ],
          },
          resolved_components: [],
        }),
      );
      renderAsk('?q=ibuprofen');

      expect(await screen.findByText('product')).toBeInTheDocument();
      expect(screen.getByText('(similar)')).toBeInTheDocument();
    });

    it('says so when the model extracted no components', async () => {
      stubApi(
        response({ interpretation: { components: [] }, resolved_components: [] }),
      );
      renderAsk('?q=anything');

      expect(await screen.findByText('(no components extracted)')).toBeInTheDocument();
    });

    it('spells out the numeric and flag filters', async () => {
      stubApi(
        response({
          interpretation: {
            components: [],
            min_yield: 70,
            max_yield: 95,
            min_conversion: 10,
            max_conversion: 90,
            reaction_smarts: '[C:1]>>[C:1]O',
            similarity_threshold: 0.8,
            use_stereochemistry: true,
            limit: 25,
          },
          resolved_components: [],
        }),
      );
      renderAsk('?q=anything');

      expect(await screen.findByText('yield ≥ 70%')).toBeInTheDocument();
      expect(screen.getByText('yield ≤ 95%')).toBeInTheDocument();
      expect(screen.getByText('conversion ≥ 10%')).toBeInTheDocument();
      expect(screen.getByText('conversion ≤ 90%')).toBeInTheDocument();
      expect(screen.getByText('reaction SMARTS [C:1]>>[C:1]O')).toBeInTheDocument();
      expect(screen.getByText('similarity threshold 0.8')).toBeInTheDocument();
      expect(screen.getByText('stereochemistry respected')).toBeInTheDocument();
      expect(screen.getByText('limit 25')).toBeInTheDocument();
    });

    // A zero bound is a real filter; only null/undefined means "not set".
    it('keeps a zero-valued filter', async () => {
      stubApi(
        response({
          interpretation: { components: [], min_yield: 0 },
          resolved_components: [],
        }),
      );
      renderAsk('?q=anything');

      expect(await screen.findByText('yield ≥ 0%')).toBeInTheDocument();
    });

    it('shows what each name resolved to, and by which resolver', async () => {
      stubApi(response());
      renderAsk('?q=reactions%20using%20benzene');

      expect(await screen.findByText('c1ccccc1')).toBeInTheDocument();
      expect(screen.getByText('via pubchem')).toBeInTheDocument();
      expect(screen.getByText('(pubchem)')).toBeInTheDocument();
    });

    // A structure the model supplied directly was never looked up, so it does
    // not belong in the resolution layer.
    it('leaves verbatim structures out of the resolution layer', async () => {
      stubApi(
        response({
          resolved_components: [
            {
              identifier: 'c1ccccc1',
              smiles: 'c1ccccc1',
              resolver: '(verbatim)',
              target: 'INPUT',
              mode: 'EXACT',
            },
          ],
        }),
      );
      renderAsk('?q=anything');

      await screen.findByText('From the model', { exact: false });
      expect(screen.queryByText(/Resolved to structures/)).not.toBeInTheDocument();
    });

    // A fresh and a cached hit from the same service are one resolver.
    it('counts a cached resolver once', async () => {
      stubApi(
        response({
          resolved_components: [
            {
              identifier: 'benzene',
              smiles: 'c1ccccc1',
              resolver: 'pubchem',
              target: 'INPUT',
              mode: 'EXACT',
            },
            {
              identifier: 'toluene',
              smiles: 'Cc1ccccc1',
              resolver: 'pubchem (cached)',
              target: 'INPUT',
              mode: 'EXACT',
            },
          ],
        }),
      );
      renderAsk('?q=anything');

      expect(await screen.findByText('(pubchem)')).toBeInTheDocument();
    });
  });

  describe('results', () => {
    it('renders the matching reactions', async () => {
      stubApi(response());
      const { container } = renderAsk('?q=reactions%20using%20benzene');

      await waitFor(() =>
        expect(container.querySelectorAll('.reaction-container')).toHaveLength(1),
      );
    });

    it('says so when nothing matched', async () => {
      stubApi(response({ results: [] }));
      renderAsk('?q=anything');

      expect(await screen.findByText(/No reactions matched/)).toBeInTheDocument();
    });

    it('reports a failed query', async () => {
      stubApi(response(), 500);
      renderAsk('?q=anything');

      expect(await screen.findByText('nl_query failed (HTTP 500)')).toBeInTheDocument();
    });
  });

  describe('dry run', () => {
    it('reads the flag from the URL and skips the search', async () => {
      stubApi(response({ dry_run: true }), 200);
      renderAsk('?q=anything&dry_run=1');

      expect(
        await screen.findByText('Dry run — search not executed'),
      ).toBeInTheDocument();
      expect(screen.getByRole('checkbox')).toBeChecked();
    });

    it('puts the flag in the URL when toggled on', async () => {
      const user = userEvent.setup();
      stubApi(response());
      renderAsk('?q=benzene');

      await user.click(screen.getByRole('checkbox'));

      await waitFor(() =>
        expect(new URLSearchParams(currentSearch).get('dry_run')).toBe('1'),
      );
      expect(new URLSearchParams(currentSearch).get('q')).toBe('benzene');
    });

    it('drops the flag from the URL when toggled off', async () => {
      const user = userEvent.setup();
      stubApi(response({ dry_run: true }));
      renderAsk('?q=benzene&dry_run=1');

      await user.click(screen.getByRole('checkbox'));

      await waitFor(() =>
        expect(new URLSearchParams(currentSearch).has('dry_run')).toBe(false),
      );
    });
  });
});
