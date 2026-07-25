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

import { fireEvent, render, screen, waitFor } from '@testing-library/react';
import { afterEach, describe, expect, it, vi } from 'vitest';
import ChartView from './ChartView';

const CHART_DATA = [
  { smiles: 'CCO', times_appearing: 12 },
  { smiles: 'CC=O', times_appearing: 4 },
];

const stubApi = (
  { stats = CHART_DATA, statsStatus = 200, svgStatus = 200 } = {} as {
    stats?: unknown;
    statsStatus?: number;
    svgStatus?: number;
  },
) => {
  const fetchMock = vi.fn(async (url: string) => {
    if (url.startsWith('/api/compound_svg')) {
      return {
        ok: svgStatus < 400,
        status: svgStatus,
        json: async () => '<svg class="mol"></svg>',
      };
    }
    return { ok: statsStatus < 400, status: statsStatus, json: async () => stats };
  });
  vi.stubGlobal('fetch', fetchMock);
  return fetchMock;
};

const renderChart = (props: Partial<React.ComponentProps<typeof ChartView>> = {}) =>
  render(
    <ChartView
      uniqueId="reactantsFrequency"
      title="Frequency of Reactants"
      apiCall="input_stats"
      role="reactant"
      datasetId="ord_dataset-1"
      {...props}
    />,
  );

const bars = (container: HTMLElement) => [...container.querySelectorAll('rect')];

afterEach(() => {
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});

describe('ChartView', () => {
  it('requests the stats for its dataset and endpoint', async () => {
    const fetchMock = stubApi();
    const { container } = renderChart();

    await waitFor(() => expect(bars(container)).toHaveLength(2));
    expect(fetchMock).toHaveBeenCalledWith(
      '/api/input_stats?dataset_id=ord_dataset-1',
      expect.objectContaining({ method: 'GET' }),
    );
  });

  it('url-encodes the dataset id', async () => {
    const fetchMock = stubApi();
    renderChart({ datasetId: 'ord dataset/1' });

    await waitFor(() =>
      expect(fetchMock).toHaveBeenCalledWith(
        '/api/input_stats?dataset_id=ord%20dataset%2F1',
        expect.anything(),
      ),
    );
  });

  it('shows the title', () => {
    stubApi();
    renderChart();
    expect(screen.getByText('Frequency of Reactants')).toBeInTheDocument();
  });

  it('draws a bar per compound', async () => {
    stubApi();
    const { container } = renderChart();
    await waitFor(() => expect(bars(container)).toHaveLength(2));
  });

  it('redraws at the collapsed size', async () => {
    stubApi();
    const { container, rerender } = renderChart({ isCollapsed: false });
    await waitFor(() =>
      expect(container.querySelector('svg')).toHaveAttribute('width', '400'),
    );

    rerender(
      <ChartView
        uniqueId="reactantsFrequency"
        title="Frequency of Reactants"
        apiCall="input_stats"
        role="reactant"
        datasetId="ord_dataset-1"
        isCollapsed
      />,
    );

    await waitFor(() =>
      expect(container.querySelector('svg')).toHaveAttribute('width', '180'),
    );
  });

  it('reports a failed stats request', async () => {
    vi.spyOn(console, 'error').mockImplementation(() => {});
    stubApi({ statsStatus: 500 });
    renderChart();

    expect(
      await screen.findByText('Failed to load chart: input_stats failed (HTTP 500)'),
    ).toBeInTheDocument();
  });

  it('draws nothing when the dataset has no stats', async () => {
    stubApi({ stats: [] });
    const { container } = renderChart();

    await waitFor(() =>
      expect(container.querySelector('.chart-view__loading')).toHaveStyle({
        visibility: 'hidden',
      }),
    );
    expect(bars(container)).toHaveLength(0);
  });

  describe('the hover tooltip', () => {
    it('shows the count and the compound drawing', async () => {
      stubApi();
      const { container } = renderChart();
      await waitFor(() => expect(bars(container)).toHaveLength(2));

      fireEvent.mouseOver(bars(container)[0], { clientX: 100, clientY: 200 });

      expect(await screen.findByText('Count: 12')).toBeInTheDocument();
      await waitFor(() =>
        expect(container.querySelector('.chart-view__svg')?.innerHTML).toBe(
          '<svg class="mol"></svg>',
        ),
      );
    });

    it('hides again on mouse out', async () => {
      stubApi();
      const { container } = renderChart();
      await waitFor(() => expect(bars(container)).toHaveLength(2));

      fireEvent.mouseOver(bars(container)[0], { clientX: 100, clientY: 200 });
      await screen.findByText('Count: 12');
      fireEvent.mouseOut(bars(container)[0]);

      expect(screen.queryByText('Count: 12')).not.toBeInTheDocument();
    });

    // The 4xx/5xx body is an HTML error page, not an SVG.
    it('does not render an error page as the drawing', async () => {
      const consoleError = vi.spyOn(console, 'error').mockImplementation(() => {});
      stubApi({ svgStatus: 500 });
      const { container } = renderChart();
      await waitFor(() => expect(bars(container)).toHaveLength(2));

      fireEvent.mouseOver(bars(container)[0], { clientX: 100, clientY: 200 });

      await waitFor(() => expect(consoleError).toHaveBeenCalled());
      expect(container.querySelector('.chart-view__svg')).toBeEmptyDOMElement();
    });
  });
});
