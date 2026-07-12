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

import React, { useCallback, useEffect, useRef, useState } from 'react';
import * as d3 from 'd3';
import { Loader, Text } from '@mantine/core';
import { ord } from 'ord-schema-protobufjs';
import { api } from '../../utils/api';
import { encodeCompound } from '../../utils/proto';
import classes from './ChartView.module.scss';

interface ChartData {
  smiles: string;
  times_appearing: number;
}

interface ChartViewProps {
  uniqueId: string;
  title: string;
  apiCall: string;
  role: string;
  datasetId: string;
}

const CHART_WIDTH = 420;
const CHART_HEIGHT = 260;

const ChartView: React.FC<ChartViewProps> = ({
  uniqueId,
  title,
  apiCall,
  role,
  datasetId,
}) => {
  const [loading, setLoading] = useState(true);
  const [fetchError, setFetchError] = useState<string | null>(null);
  const [inputsData, setInputsData] = useState<ChartData[]>([]);
  const [showTooltip, setShowTooltip] = useState<'visible' | 'hidden'>('hidden');
  const [currentTimesAppearing, setCurrentTimesAppearing] = useState(0);
  const [tooltipOffsetHorizontal, setTooltipOffsetHorizontal] = useState(0);
  const [tooltipOffsetVertical, setTooltipOffsetVertical] = useState(0);
  const [showSmiles, setShowSmiles] = useState(false);
  const [molHtml, setMolHtml] = useState<string | null>(null);
  const [molLoading, setMolLoading] = useState(true);

  const svgRef = useRef<SVGSVGElement>(null);

  const getMolHtml = useCallback(async (smiles: string): Promise<string> => {
    const binary = encodeCompound({
      identifiers: [
        {
          type: ord.CompoundIdentifier.CompoundIdentifierType.SMILES,
          value: smiles,
        },
      ],
    });

    const response = await api.post<string>('/compound_svg', binary, {
      headers: { 'Content-Type': 'application/x-protobuf' },
    });
    return response.data;
  }, []);

  const createChart = useCallback(
    (data: ChartData[], width: number, height: number) => {
      if (!svgRef.current) return;

      // Clear previous chart
      d3.select(svgRef.current).selectAll('*').remove();

      const marginTop = 20;
      const marginRight = 20;
      const marginBottom = 30;
      const marginLeft = 40;

      // X-axis
      const x = d3
        .scaleBand()
        .domain(
          d3.groupSort(
            data,
            ([d]) => -d.times_appearing,
            d => d.smiles,
          ),
        )
        .range([marginLeft, width - marginRight])
        .padding(0.1);

      // Y-axis
      const y = d3
        .scaleLinear()
        .domain([0, d3.max(data, d => d.times_appearing) || 0])
        .range([height - marginBottom, marginTop]);

      // SVG setup
      const svg = d3
        .select(svgRef.current)
        .attr('width', width)
        .attr('height', height)
        .attr('viewBox', [0, 0, width, height])
        .attr('style', 'max-width: 100%; height: auto;');

      // Create bars in the shared primary blue.
      svg
        .append('g')
        .attr('fill', '#3c78d8')
        .selectAll('rect')
        .data(data)
        .join('rect')
        .attr('x', d => x(d.smiles) || 0)
        .attr('y', d => y(d.times_appearing))
        .attr('height', d => y(0) - y(d.times_appearing))
        .attr('width', x.bandwidth())
        .attr('style', 'cursor: pointer;')
        .on('mouseover', (event: MouseEvent, d: ChartData) => {
          setCurrentTimesAppearing(d.times_appearing);
          setTooltipOffsetHorizontal(event.clientX);
          setTooltipOffsetVertical(
            role === 'product' ? event.clientY - 240 : event.clientY - 140,
          );
          setShowTooltip('visible');
          setShowSmiles(true);
          setMolLoading(true);

          getMolHtml(d.smiles)
            .then(result => {
              setMolLoading(false);
              setMolHtml(result);
            })
            .catch(error => {
              console.error('Error fetching molecule SVG:', error);
              setMolLoading(false);
              setMolHtml(null);
            });
        })
        .on('mouseout', () => {
          setShowTooltip('hidden');
          setMolHtml(null);
          setShowSmiles(false);
          setMolLoading(true);
        });

      // X-axis formatting
      svg
        .append('g')
        .attr('transform', `translate(0,${height - marginBottom})`)
        .call(
          d3
            .axisBottom(x)
            .ticks(width / 80)
            .tickSizeOuter(0)
            .tickFormat(() => ''),
        )
        .call(g =>
          g
            .append('text')
            .attr('x', width)
            .attr('y', marginBottom - 4)
            .attr('fill', 'currentColor')
            .attr('text-anchor', 'end')
            .text('Molecules (Hover to view) →'),
        );

      // Y-axis formatting
      svg
        .append('g')
        .attr('transform', `translate(${marginLeft},0)`)
        .call(d3.axisLeft(y).tickFormat(y => (y as number).toFixed()))
        .call(g => g.select('.domain').remove())
        .call(g =>
          g
            .append('text')
            .attr('x', -marginLeft)
            .attr('y', 10)
            .attr('fill', 'currentColor')
            .attr('text-anchor', 'start')
            .text('↑ Frequency (no. of occurrences)'),
        );
    },
    [getMolHtml, role],
  );

  // Fetch chart data when the endpoint or dataset changes. The AbortController
  // guards against a stale response from the previous datasetId racing in
  // after the new fetch has started.
  useEffect(() => {
    const controller = new AbortController();
    setLoading(true);
    setFetchError(null);
    api
      .get<ChartData[]>(`/${apiCall}`, {
        params: { dataset_id: datasetId },
        signal: controller.signal,
      })
      .then(response => {
        setLoading(false);
        setInputsData(response.data);
      })
      .catch((error: Error) => {
        if (error.name === 'CanceledError' || error.name === 'AbortError') return;
        console.error(`Error fetching ${apiCall}:`, error);
        setLoading(false);
        setFetchError(error.message);
      });
    return () => controller.abort();
  }, [apiCall, datasetId]);

  useEffect(() => {
    if (inputsData.length > 0) {
      createChart(inputsData, CHART_WIDTH, CHART_HEIGHT);
    }
  }, [inputsData, createChart]);

  return (
    <div className={classes.chartView}>
      <div className={classes.titleAndChart}>
        <Text
          className={classes.title}
          fw={500}
        >
          {title}
        </Text>

        <svg
          ref={svgRef}
          id={uniqueId}
          style={{ visibility: loading || fetchError ? 'hidden' : 'visible' }}
        />
      </div>

      {fetchError ? (
        <div className={classes.error}>Failed to load chart: {fetchError}</div>
      ) : (
        <div
          className={classes.loading}
          style={{ visibility: loading ? 'visible' : 'hidden' }}
        >
          <Loader size="sm" />
        </div>
      )}

      {showSmiles && (
        <div
          className={classes.tooltip}
          style={{
            top: `${tooltipOffsetVertical}px`,
            left: `${tooltipOffsetHorizontal}px`,
            visibility: showTooltip,
          }}
        >
          <pre>Count: {currentTimesAppearing}</pre>
          <div
            className={classes.molSvg}
            dangerouslySetInnerHTML={{ __html: molHtml || '' }}
          />
          <div
            className={classes.molLoading}
            style={{ visibility: molLoading ? 'visible' : 'hidden' }}
          >
            <Loader size="sm" />
          </div>
        </div>
      )}
    </div>
  );
};

export default ChartView;
