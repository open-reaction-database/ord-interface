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
import { MemoryRouter, Route, Routes, useLocation } from 'react-router-dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import SearchResults from './SearchResults';
import type { ReactionData, SearchResult } from '../../types/search';

const SEARCH = '?dataset_id=ord_dataset-1';

const results = (...ids: string[]): SearchResult[] =>
  ids.map(id => ({
    reaction_id: id,
    dataset_id: 'ord_dataset-1',
    proto: '',
    data: {} as ReactionData,
  }));

let selectedSetSearch = '';

const SelectedSetPage = () => {
  selectedSetSearch = useLocation().search;
  return <div>selected set page</div>;
};

const renderResults = (searchResults: SearchResult[], search = SEARCH) =>
  render(
    <MemoryRouter initialEntries={[`/search${search}`]}>
      <Routes>
        <Route
          path="/search"
          element={<SearchResults searchResults={searchResults} />}
        />
        <Route
          path="/selected-set"
          element={<SelectedSetPage />}
        />
      </Routes>
    </MemoryRouter>,
  );

const checkboxes = () => screen.getAllByLabelText('Select reaction');

beforeEach(() => {
  selectedSetSearch = '';
  vi.stubGlobal(
    'fetch',
    vi.fn(async () => ({ ok: true, status: 200, text: async () => '' })),
  );
});

afterEach(() => {
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});

describe('SearchResults', () => {
  it('renders no table without results', () => {
    renderResults([]);
    expect(screen.queryByText('Search Results')).not.toBeInTheDocument();
    expect(screen.queryByLabelText('Select reaction')).not.toBeInTheDocument();
  });

  it('renders a card per result', async () => {
    renderResults(results('ord-1', 'ord-2'));
    expect(await screen.findAllByLabelText('Select reaction')).toHaveLength(2);
  });

  it('offers a shareable link and a download', async () => {
    renderResults(results('ord-1'));
    expect(await screen.findByText('Shareable Link')).toBeInTheDocument();
    expect(
      screen.getByRole('button', { name: 'Download All Search Results' }),
    ).toBeEnabled();
  });

  it('opens the download modal for every result, not just the selected ones', async () => {
    const user = userEvent.setup();
    renderResults(results('ord-1', 'ord-2'));

    await user.click(
      await screen.findByRole('button', { name: 'Download All Search Results' }),
    );

    expect(screen.getByText('Download Results')).toBeInTheDocument();
  });

  describe('selection', () => {
    it('counts the selected reactions', async () => {
      const user = userEvent.setup();
      renderResults(results('ord-1', 'ord-2'));

      await user.click((await screen.findAllByLabelText('Select reaction'))[0]);
      expect(screen.getByText('View 1 selected reactions')).toBeInTheDocument();

      await user.click(checkboxes()[1]);
      expect(screen.getByText('View 2 selected reactions')).toBeInTheDocument();
    });

    it('drops a deselected reaction', async () => {
      const user = userEvent.setup();
      renderResults(results('ord-1', 'ord-2'));

      const boxes = await screen.findAllByLabelText('Select reaction');
      await user.click(boxes[0]);
      await user.click(boxes[0]);

      expect(screen.queryByText(/selected reactions/)).not.toBeInTheDocument();
    });

    it('carries the selection to the selected-set page', async () => {
      const user = userEvent.setup();
      renderResults(results('ord-1', 'ord-2'));

      await user.click((await screen.findAllByLabelText('Select reaction'))[1]);
      await user.click(screen.getByText('View 1 selected reactions'));

      expect(screen.getByText('selected set page')).toBeInTheDocument();
      expect(new URLSearchParams(selectedSetSearch).getAll('reaction_id')).toEqual([
        'ord-2',
      ]);
    });

    // The selection has to survive navigating to the selected-set page and back.
    it('stores the selection alongside the query that produced it', async () => {
      const user = userEvent.setup();
      renderResults(results('ord-1'));

      await user.click((await screen.findAllByLabelText('Select reaction'))[0]);
      await user.click(screen.getByText('View 1 selected reactions'));

      expect(JSON.parse(localStorage.getItem('storedSet')!)).toEqual({
        query: SEARCH,
        reactions: ['ord-1'],
      });
    });

    it('restores a stored selection for the same query', async () => {
      localStorage.setItem(
        'storedSet',
        JSON.stringify({ query: SEARCH, reactions: ['ord-2'] }),
      );
      renderResults(results('ord-1', 'ord-2'));

      const boxes = await screen.findAllByLabelText('Select reaction');
      expect(boxes[0]).not.toBeChecked();
      expect(boxes[1]).toBeChecked();
    });

    it('ignores a selection stored under a different query', async () => {
      localStorage.setItem(
        'storedSet',
        JSON.stringify({ query: '?dataset_id=other', reactions: ['ord-1'] }),
      );
      renderResults(results('ord-1'));

      expect((await screen.findAllByLabelText('Select reaction'))[0]).not.toBeChecked();
    });

    it('survives a corrupt stored selection', async () => {
      const consoleError = vi.spyOn(console, 'error').mockImplementation(() => {});
      localStorage.setItem('storedSet', 'not json');
      renderResults(results('ord-1'));

      expect(await screen.findAllByLabelText('Select reaction')).toHaveLength(1);
      expect(consoleError).toHaveBeenCalled();
    });
  });
});
