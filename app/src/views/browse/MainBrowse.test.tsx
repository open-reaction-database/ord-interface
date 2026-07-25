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
import { afterEach, describe, expect, it, vi } from 'vitest';
import MainBrowse from './MainBrowse';

const dataset = (overrides = {}) => ({
  dataset_id: 'ord_dataset-1',
  name: 'Suzuki couplings',
  description: 'Reactions from https://doi.org/10.1021/jacs.8b01523',
  num_reactions: 1234,
  submitted_at: '2021-05-04',
  ...overrides,
});

const stubDatasets = (datasets: unknown[]) => {
  const fetchMock = vi.fn(async () => ({ json: async () => datasets }));
  vi.stubGlobal('fetch', fetchMock);
  return fetchMock;
};

const renderBrowse = () =>
  render(
    <MemoryRouter>
      <MainBrowse />
    </MemoryRouter>,
  );

afterEach(() => {
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});

describe('MainBrowse', () => {
  it('requests the dataset list', async () => {
    const fetchMock = stubDatasets([dataset()]);
    renderBrowse();

    expect(await screen.findByText('ord_dataset-1')).toBeInTheDocument();
    expect(fetchMock).toHaveBeenCalledWith('/api/datasets', { method: 'GET' });
  });

  it('shows a spinner until the datasets arrive', () => {
    stubDatasets([dataset()]);
    const { container } = renderBrowse();
    expect(container.querySelector('.spinner-main')).toBeInTheDocument();
  });

  it('keeps the spinner up when the list comes back empty', async () => {
    stubDatasets([]);
    const { container } = renderBrowse();
    expect(container.querySelector('.spinner-main')).toBeInTheDocument();
  });

  it('links each dataset to its detail page', async () => {
    stubDatasets([dataset()]);
    renderBrowse();

    expect(await screen.findByRole('link', { name: 'ord_dataset-1' })).toHaveAttribute(
      'href',
      '/dataset/ord_dataset-1',
    );
  });

  it('groups the reaction count', async () => {
    stubDatasets([dataset({ num_reactions: 1234567 })]);
    renderBrowse();

    expect(await screen.findByText('1,234,567')).toBeInTheDocument();
  });

  // Names and descriptions are free text that often carries a DOI or URL.
  it('links references in the name and description', async () => {
    stubDatasets([dataset({ name: 'From 10.1021/jacs.8b01523' })]);
    renderBrowse();

    const links = await screen.findAllByRole('link', {
      name: /10\.1021\/jacs\.8b01523/,
    });
    expect(links[0]).toHaveAttribute('href', 'https://doi.org/10.1021/jacs.8b01523');
  });

  describe('the details modal', () => {
    it('opens with the submission date and reaction count', async () => {
      const user = userEvent.setup();
      stubDatasets([dataset()]);
      renderBrowse();

      await user.click(
        await screen.findByRole('button', { name: 'Show dataset details' }),
      );

      expect(screen.getByText('Submitted:')).toBeInTheDocument();
      expect(screen.getByText('2021-05-04')).toBeInTheDocument();
      expect(screen.getByText('Reactions:')).toBeInTheDocument();
      expect(screen.getAllByText('1,234')).not.toHaveLength(0);
    });

    it('shows a dash for an unknown submission date', async () => {
      const user = userEvent.setup();
      stubDatasets([dataset({ submitted_at: null })]);
      renderBrowse();

      await user.click(
        await screen.findByRole('button', { name: 'Show dataset details' }),
      );

      expect(screen.getByText('—')).toBeInTheDocument();
    });

    it('closes again', async () => {
      const user = userEvent.setup();
      stubDatasets([dataset()]);
      renderBrowse();

      await user.click(
        await screen.findByRole('button', { name: 'Show dataset details' }),
      );
      await user.click(screen.getByText('✕'));

      expect(screen.queryByText('Submitted:')).not.toBeInTheDocument();
    });
  });

  it('survives a failed request', async () => {
    const consoleError = vi.spyOn(console, 'error').mockImplementation(() => {});
    vi.stubGlobal(
      'fetch',
      vi.fn(async () => Promise.reject(new Error('offline'))),
    );

    const { container } = renderBrowse();

    await vi.waitFor(() => expect(consoleError).toHaveBeenCalled());
    expect(container.querySelector('.spinner-main')).toBeInTheDocument();
  });
});
