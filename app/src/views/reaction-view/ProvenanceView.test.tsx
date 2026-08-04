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
import ProvenanceView from './ProvenanceView';
import type { ReactionProvenanceData } from '../../types/search';

const renderProvenance = (provenance: unknown) =>
  render(<ProvenanceView provenance={provenance as ReactionProvenanceData} />);

const fields = (container: HTMLElement, selector: string): Record<string, string> => {
  const details = container.querySelector(selector);
  if (!details) return {};
  const labels = [...details.querySelectorAll('.label')].map(
    el => el.textContent ?? '',
  );
  const values = [...details.querySelectorAll('.value')].map(
    el => el.textContent ?? '',
  );
  return Object.fromEntries(labels.map((label, index) => [label, values[index] ?? '']));
};

describe('ProvenanceView', () => {
  it('renders an empty view without provenance', () => {
    const { container } = renderProvenance(undefined);
    expect(container.querySelector('.provenance-view')).toBeInTheDocument();
    expect(container.querySelector('.experimenter')).toBeNull();
  });

  it('shows the experimenter when one was recorded', () => {
    const { container } = renderProvenance({
      experimenter: {
        username: 'alovelace',
        name: 'Ada Lovelace',
        orcid: '0000-0001-2345-6789',
        organization: 'ORD',
        email: 'ada@example.com',
      },
    });

    expect(fields(container, '.experimenter')).toEqual({
      Username: 'alovelace',
      Name: 'Ada Lovelace',
      ORCID: '0000-0001-2345-6789',
      Organization: 'ORD',
      Email: 'ada@example.com',
    });
  });

  it('omits the experimenter section when none was recorded', () => {
    renderProvenance({ city: 'Cambridge' });
    expect(screen.queryByText('Experimenter')).not.toBeInTheDocument();
  });

  it('shows the city, DOI and patent', () => {
    const { container } = renderProvenance({
      city: 'Cambridge',
      doi: '10.1021/jacs.8b01523',
      patent: 'US1234567',
    });

    expect(fields(container, '.details')).toEqual({
      City: 'Cambridge',
      DOI: '10.1021/jacs.8b01523',
      Patent: 'US1234567',
    });
  });

  it('links the publication URL', () => {
    renderProvenance({ publicationUrl: 'https://example.com/paper' });
    expect(
      screen.getByRole('link', { name: 'https://example.com/paper' }),
    ).toHaveAttribute('href', 'https://example.com/paper');
  });

  it('omits the experiment start when the timestamp is empty', () => {
    renderProvenance({ experimentStart: { value: '' }, city: 'Cambridge' });
    expect(screen.queryByText('Experiment Start')).not.toBeInTheDocument();
  });

  it('shows the experiment start as a date', () => {
    renderProvenance({ experimentStart: { value: '2021-05-04T12:00:00Z' } });
    expect(screen.getByText('Experiment Start')).toBeInTheDocument();
  });
});
