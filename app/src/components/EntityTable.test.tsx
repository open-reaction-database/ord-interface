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
import { describe, expect, it } from 'vitest';
import EntityTable from './EntityTable';

interface Row {
  id: string;
  name: string;
  organization?: string | null;
}

const rows = (count: number): Row[] =>
  Array.from({ length: count }, (_, i) => ({
    id: `row-${i}`,
    name: `Name ${i}`,
  }));

const renderTable = (props: Partial<React.ComponentProps<typeof EntityTable<Row>>>) =>
  render(
    <EntityTable<Row>
      tableData={[]}
      {...props}
    >
      {entities => (
        <ul>
          {entities.map(entity => (
            <li key={entity.id}>{entity.id}</li>
          ))}
        </ul>
      )}
    </EntityTable>,
  );

// The ids the render prop was handed on the current page.
const visibleIds = (): string[] =>
  screen.queryAllByRole('listitem').map(item => item.textContent ?? '');

describe('EntityTable', () => {
  it('renders nothing without data', () => {
    const { container } = renderTable({ tableData: [] });
    expect(container).toBeEmptyDOMElement();
  });

  it('shows the title', () => {
    renderTable({ tableData: rows(1), title: 'Search Results' });
    expect(screen.getByText('Search Results')).toBeInTheDocument();
  });

  it('shows the search box by default and hides it on request', () => {
    const { unmount } = renderTable({ tableData: rows(1) });
    expect(screen.getByLabelText('Search:')).toBeInTheDocument();
    unmount();

    renderTable({ tableData: rows(1), displaySearch: false });
    expect(screen.queryByLabelText('Search:')).not.toBeInTheDocument();
  });

  it('hands the render prop only the first page', () => {
    renderTable({ tableData: rows(25) });
    expect(visibleIds()).toEqual(rows(10).map(row => row.id));
    expect(screen.getByText(/of 25 entries/)).toBeInTheDocument();
  });

  it('paginates forward and back', async () => {
    const user = userEvent.setup();
    renderTable({ tableData: rows(25) });

    await user.click(screen.getByText('Next'));
    expect(visibleIds()[0]).toBe('row-10');

    await user.click(screen.getByText('Last'));
    expect(visibleIds()).toEqual(['row-20', 'row-21', 'row-22', 'row-23', 'row-24']);

    await user.click(screen.getByText('First'));
    expect(visibleIds()[0]).toBe('row-0');
  });

  // "Next" and "Previous" wrap around; the numbered links next to the current
  // page do not.
  it('wraps around at the ends', async () => {
    const user = userEvent.setup();
    renderTable({ tableData: rows(25) });

    await user.click(screen.getByText('Previous'));
    expect(visibleIds()[0]).toBe('row-20');

    await user.click(screen.getByText('Next'));
    expect(visibleIds()[0]).toBe('row-0');
  });

  it('resizes the page and returns to the first page', async () => {
    const user = userEvent.setup();
    renderTable({ tableData: rows(25) });

    await user.click(screen.getByText('Last'));
    expect(visibleIds()[0]).toBe('row-20');

    await user.selectOptions(screen.getByRole('combobox'), '25');
    expect(visibleIds()).toHaveLength(25);
  });

  // Filtering from a later page used to leave the user on a page number past
  // the end of the narrowed results, showing an empty list.
  it('returns to the first page when the filter changes', async () => {
    const user = userEvent.setup();
    renderTable({ tableData: rows(25) });

    await user.click(screen.getByText('Last'));
    expect(visibleIds()[0]).toBe('row-20');

    await user.type(screen.getByLabelText('Search:'), 'row-1');
    expect(visibleIds()).toEqual([
      'row-1',
      'row-10',
      'row-11',
      'row-12',
      'row-13',
      'row-14',
      'row-15',
      'row-16',
      'row-17',
      'row-18',
    ]);
  });

  it('returns to the first page when the filter is cleared', async () => {
    const user = userEvent.setup();
    renderTable({ tableData: rows(25) });

    await user.type(screen.getByLabelText('Search:'), 'row-1');
    await user.click(screen.getByText('Last'));
    expect(visibleIds()).toEqual(['row-19']);

    await user.clear(screen.getByLabelText('Search:'));
    expect(visibleIds()[0]).toBe('row-0');
  });

  it('filters on any field, case-insensitively', async () => {
    const user = userEvent.setup();
    renderTable({
      tableData: [
        { id: 'a', name: 'Suzuki', organization: 'MIT' },
        { id: 'b', name: 'Heck', organization: 'Stanford' },
      ],
    });

    await user.type(screen.getByLabelText('Search:'), 'suzuki');
    expect(visibleIds()).toEqual(['a']);

    await user.clear(screen.getByLabelText('Search:'));
    await user.type(screen.getByLabelText('Search:'), 'stanford');
    expect(visibleIds()).toEqual(['b']);
  });

  // Terms are ANDed: every whitespace-separated term must match some field.
  it('requires every search term to match', async () => {
    const user = userEvent.setup();
    renderTable({
      tableData: [
        { id: 'a', name: 'Suzuki', organization: 'MIT' },
        { id: 'b', name: 'Suzuki', organization: 'Stanford' },
      ],
    });

    await user.type(screen.getByLabelText('Search:'), 'suzuki mit');
    expect(visibleIds()).toEqual(['a']);
  });

  it('ignores extra whitespace between terms', async () => {
    const user = userEvent.setup();
    renderTable({ tableData: rows(3) });

    await user.type(screen.getByLabelText('Search:'), '  row-1   ');
    expect(visibleIds()).toEqual(['row-1']);
  });

  it('reports an empty result set rather than falling back to everything', async () => {
    const user = userEvent.setup();
    renderTable({ tableData: rows(3) });

    await user.type(screen.getByLabelText('Search:'), 'no-such-row');
    expect(visibleIds()).toEqual([]);
    expect(screen.getByText(/of 0 entries/)).toBeInTheDocument();
  });

  it('matches an absent field as the literal "null"', async () => {
    const user = userEvent.setup();
    renderTable({
      tableData: [
        { id: 'a', name: 'Suzuki', organization: null },
        { id: 'b', name: 'Heck', organization: 'Stanford' },
      ],
    });

    await user.type(screen.getByLabelText('Search:'), 'null');
    expect(visibleIds()).toEqual(['a']);
  });

  it('follows a change of data', () => {
    const { rerender } = render(
      <EntityTable<Row> tableData={rows(2)}>
        {entities => (
          <ul>
            {entities.map(entity => (
              <li key={entity.id}>{entity.id}</li>
            ))}
          </ul>
        )}
      </EntityTable>,
    );
    expect(visibleIds()).toHaveLength(2);

    rerender(
      <EntityTable<Row> tableData={rows(5)}>
        {entities => (
          <ul>
            {entities.map(entity => (
              <li key={entity.id}>{entity.id}</li>
            ))}
          </ul>
        )}
      </EntityTable>,
    );
    expect(visibleIds()).toHaveLength(5);
  });
});
