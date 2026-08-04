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
import { MemoryRouter, Route, Routes } from 'react-router-dom';
import { afterEach, describe, expect, it, vi } from 'vitest';
import MainDatasetView from './MainDatasetView';

const encodedReaction = (reactionId: string): string => {
  const reaction = new reaction_pb.Reaction();
  reaction.setReactionId(reactionId);
  return btoa(String.fromCharCode(...reaction.serializeBinary()));
};

interface ApiOverrides {
  dataset?: { status?: number; body?: unknown };
  reactions?: string[];
  resultStatus?: number;
}

const stubApi = ({
  dataset,
  reactions = ['ord-1'],
  resultStatus = 200,
}: ApiOverrides) => {
  const fetchMock = vi.fn(async (url: string) => {
    if (url.startsWith('/api/dataset?')) {
      const status = dataset?.status ?? 200;
      return {
        ok: status < 400,
        status,
        json: async () =>
          dataset?.body ?? {
            dataset_id: 'ord_dataset-1',
            name: 'Suzuki couplings',
            description: 'from the literature',
            num_reactions: 500,
          },
      };
    }
    if (url.startsWith('/api/submit_query')) {
      return { ok: true, status: 200, json: async () => 'task-1' };
    }
    if (url.startsWith('/api/fetch_query_result')) {
      return {
        ok: resultStatus < 400,
        status: resultStatus,
        json: async () =>
          reactions.map(id => ({
            reaction_id: id,
            dataset_id: 'ord_dataset-1',
            proto: encodedReaction(id),
          })),
      };
    }
    // Chart data and per-card summaries.
    if (url.startsWith('/api/reaction_summary')) {
      return { ok: true, status: 200, text: async () => '' };
    }
    return { ok: true, status: 200, json: async () => [] };
  });
  vi.stubGlobal('fetch', fetchMock);
  return fetchMock;
};

const renderDataset = (datasetId = 'ord_dataset-1') => {
  const client = new QueryClient({
    defaultOptions: { queries: { retry: false, gcTime: 0 } },
  });
  return render(
    <QueryClientProvider client={client}>
      <MemoryRouter initialEntries={[`/dataset/${datasetId}`]}>
        <Routes>
          <Route
            path="/dataset/:datasetId"
            element={<MainDatasetView />}
          />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
};

afterEach(() => {
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});

describe('MainDatasetView', () => {
  it('shows the dataset metadata', async () => {
    stubApi({});
    renderDataset();

    expect(await screen.findByText('Suzuki couplings')).toBeInTheDocument();
    expect(screen.getByText('from the literature')).toBeInTheDocument();
    expect(screen.getByText('500')).toBeInTheDocument();
  });

  it('fills in placeholders for missing metadata', async () => {
    stubApi({ dataset: { body: { dataset_id: 'ord_dataset-1' } } });
    renderDataset();

    expect(await screen.findByText('(no name)')).toBeInTheDocument();
    expect(screen.getByText('(no description)')).toBeInTheDocument();
  });

  it('reports failed metadata without hiding the rest of the page', async () => {
    stubApi({ dataset: { status: 500 } });
    renderDataset();

    expect(
      await screen.findByText(
        'Failed to load dataset metadata: dataset metadata failed (HTTP 500)',
      ),
    ).toBeInTheDocument();
  });

  it('searches the dataset for its reactions', async () => {
    const fetchMock = stubApi({});
    renderDataset();

    await waitFor(() =>
      expect(fetchMock).toHaveBeenCalledWith(
        '/api/submit_query?dataset_id=ord_dataset-1&limit=100',
        undefined,
      ),
    );
  });

  it('lists the reactions it found', async () => {
    stubApi({ reactions: ['ord-1', 'ord-2'] });
    const { container } = renderDataset();

    await waitFor(() =>
      expect(container.querySelectorAll('.reaction-container')).toHaveLength(2),
    );
  });

  // The search is capped at 100, so a bigger dataset is shown as a sample.
  it('marks the results as a sample when the dataset is larger', async () => {
    stubApi({ reactions: ['ord-1'] });
    renderDataset();

    expect(
      await screen.findByText('100 Reactions From This Dataset (Sample)'),
    ).toBeInTheDocument();
  });

  it('says so when the dataset has no reactions', async () => {
    stubApi({ reactions: [] });
    renderDataset();

    expect(
      await screen.findByText('This dataset contains no reactions.'),
    ).toBeInTheDocument();
  });

  it('reports a failed reaction search', async () => {
    stubApi({ resultStatus: 500 });
    renderDataset();

    expect(await screen.findByText(/Failed to load reactions:/)).toBeInTheDocument();
  });

  it('toggles the chart panel between collapsed and expanded', async () => {
    const user = userEvent.setup();
    stubApi({});
    renderDataset();

    expect(await screen.findByTitle('Expand')).toBeInTheDocument();
    await user.click(screen.getByRole('button', { name: /keyboard_double_arrow/ }));
    expect(screen.getByTitle('Collapse')).toBeInTheDocument();
  });

  it('renders a chart for reactants and one for products', async () => {
    stubApi({});
    renderDataset();

    expect(await screen.findByText('Frequency of Reactants')).toBeInTheDocument();
    expect(screen.getByText('Frequency of Products')).toBeInTheDocument();
  });
});
