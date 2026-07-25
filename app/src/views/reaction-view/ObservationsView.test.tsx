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
import { describe, expect, it } from 'vitest';
import ObservationsView from './ObservationsView';
import type { ReactionObservationData } from '../../types/search';

const renderObservations = (observations: unknown[]) =>
  render(<ObservationsView observations={observations as ReactionObservationData[]} />);

// The rows are a flat grid of .value cells, two per observation.
const rows = (container: HTMLElement): string[][] => {
  const cells = [...container.querySelectorAll('.value')].map(
    el => el.textContent ?? '',
  );
  return cells.reduce<string[][]>((pairs, cell, index) => {
    if (index % 2 === 0) pairs.push([cell]);
    else pairs[pairs.length - 1].push(cell);
    return pairs;
  }, []);
};

describe('ObservationsView', () => {
  it('renders nothing without observations', () => {
    const { container } = renderObservations([]);
    expect(container).toBeEmptyDOMElement();
  });

  it('heads the table with Time and Comment', () => {
    renderObservations([{ comment: 'turned yellow' }]);
    expect(screen.getByText('Time')).toBeInTheDocument();
    expect(screen.getByText('Comment')).toBeInTheDocument();
  });

  it('renders a row per observation, in order', () => {
    const { container } = renderObservations([
      { time: { value: 5, units: 2 }, comment: 'turned yellow' },
      { time: { value: 2, units: 1 }, comment: 'precipitate formed' },
    ]);

    expect(rows(container)).toEqual([
      ['5 minute(s)', 'turned yellow'],
      ['2 hour(s)', 'precipitate formed'],
    ]);
  });

  it('leaves the time blank when none was recorded', () => {
    const { container } = renderObservations([{ comment: 'turned yellow' }]);
    expect(rows(container)).toEqual([['', 'turned yellow']]);
  });
});
