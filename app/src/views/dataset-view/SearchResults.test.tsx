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
import { MemoryRouter } from 'react-router-dom';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import SearchResults from './SearchResults';
import type { ReactionData, SearchResult } from '../../types/search';

const results = (count: number): SearchResult[] =>
  Array.from({ length: count }, (_, index) => ({
    reaction_id: `ord-${index}`,
    dataset_id: 'ord_dataset-1',
    proto: '',
    data: {} as ReactionData,
  }));

const renderResults = (searchResults: SearchResult[], isOverflow = false) =>
  render(
    <MemoryRouter>
      <SearchResults
        searchResults={searchResults}
        isOverflow={isOverflow}
      />
    </MemoryRouter>,
  );

beforeEach(() => {
  vi.stubGlobal(
    'fetch',
    vi.fn(async () => ({ ok: true, status: 200, text: async () => '' })),
  );
});

afterEach(() => {
  vi.unstubAllGlobals();
});

describe('dataset SearchResults', () => {
  it('renders an empty shell without results', () => {
    const { container } = renderResults([]);
    expect(container.querySelector('.search-results')).toBeEmptyDOMElement();
  });

  it('counts the reactions in the title', async () => {
    renderResults(results(3));
    expect(
      await screen.findByText('Reactions in this Dataset (3 Reactions)'),
    ).toBeInTheDocument();
  });

  // A dataset larger than the 100-result cap is shown as a sample, and the
  // title has to say so rather than implying the dataset is that small.
  it('says so when the results are only a sample', async () => {
    renderResults(results(3), true);
    expect(
      await screen.findByText('100 Reactions From This Dataset (Sample)'),
    ).toBeInTheDocument();
  });

  // Selection belongs to search, not to browsing a dataset.
  it('renders the cards without selection checkboxes', async () => {
    renderResults(results(2));
    await screen.findByText(/Reactions in this Dataset/);
    expect(screen.queryByLabelText('Select reaction')).not.toBeInTheDocument();
  });

  it('opens the download modal', async () => {
    const user = userEvent.setup();
    renderResults(results(2));

    await user.click(
      await screen.findByRole('button', { name: 'Download All Search Results' }),
    );

    expect(screen.getByText('Download Results')).toBeInTheDocument();
  });
});
