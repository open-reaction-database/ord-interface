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

import { render } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import type { RecordEvent } from 'ord-schema/proto/reaction_pb';
import { describe, expect, it } from 'vitest';
import EventsView from './EventsView';

const renderEvents = (events: unknown[]) =>
  render(<EventsView events={events as RecordEvent.AsObject[]} />);

const fields = (container: HTMLElement): Record<string, string> => {
  const details = container.querySelector('.details');
  if (!details) return {};
  const labels = [...details.querySelectorAll('.label')].map(
    el => el.textContent ?? '',
  );
  const values = [...details.querySelectorAll('.value')].map(
    el => el.textContent ?? '',
  );
  return Object.fromEntries(labels.map((label, index) => [label, values[index] ?? '']));
};

const created = (details: string, person: Record<string, string> = {}) => ({
  time: { value: '2021-05-04T12:00:00Z' },
  details,
  person,
});

describe('EventsView', () => {
  it('renders nothing without events', () => {
    const { container } = renderEvents([]);
    expect(container).toBeEmptyDOMElement();
  });

  it('shows the details and every recorded person field', () => {
    const { container } = renderEvents([
      created('record created', {
        username: 'alovelace',
        name: 'Ada Lovelace',
        orcid: '0000-0001-2345-6789',
        organization: 'ORD',
        email: 'ada@example.com',
      }),
    ]);

    expect(fields(container)).toEqual({
      Details: 'record created',
      Username: 'alovelace',
      Name: 'Ada Lovelace',
      ORCID: '0000-0001-2345-6789',
      Organization: 'ORD',
      Email: 'ada@example.com',
    });
  });

  it('omits the person fields that were not recorded', () => {
    const { container } = renderEvents([created('record created', { name: 'Ada' })]);
    expect(fields(container)).toEqual({ Details: 'record created', Name: 'Ada' });
  });

  it('opens on the first event and switches on click', async () => {
    const user = userEvent.setup();
    const { container } = renderEvents([created('first edit'), created('second edit')]);
    expect(fields(container)).toEqual({ Details: 'first edit' });

    const tabs = [...container.querySelectorAll('.tab')];
    await user.click(tabs[1]);

    expect(fields(container)).toEqual({ Details: 'second edit' });
    expect(tabs[1]).toHaveClass('selected');
  });

  it('leaves the tab unlabeled when the event has no timestamp', () => {
    const { container } = renderEvents([{ details: 'undated', person: {} }]);
    expect(container.querySelector('.tab')).toBeEmptyDOMElement();
  });
});
