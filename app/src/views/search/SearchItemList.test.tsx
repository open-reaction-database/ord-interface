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
import { describe, expect, it, vi } from 'vitest';
import SearchItemList from './SearchItemList';

const renderList = (itemList: string[] = []) => {
  const onUpdateItemList = vi.fn();
  render(
    <SearchItemList
      title="Reaction IDs"
      itemList={itemList}
      onUpdateItemList={onUpdateItemList}
    />,
  );
  return { onUpdateItemList, input: screen.getByRole('textbox') };
};

describe('SearchItemList', () => {
  it('shows the title and the items it was given', () => {
    renderList(['ord-1', 'ord-2']);
    expect(screen.getByText('Reaction IDs')).toBeInTheDocument();
    expect(screen.getByText('ord-1')).toBeInTheDocument();
    expect(screen.getByText('ord-2')).toBeInTheDocument();
  });

  it('adds an item and reports the new list', async () => {
    const user = userEvent.setup();
    const { onUpdateItemList, input } = renderList(['ord-1']);

    await user.type(input, 'ord-2');
    await user.click(screen.getByRole('button', { name: /add/i }));

    expect(onUpdateItemList).toHaveBeenCalledWith(['ord-1', 'ord-2']);
    expect(screen.getByText('ord-2')).toBeInTheDocument();
  });

  it('clears the input after adding', async () => {
    const user = userEvent.setup();
    const { input } = renderList();

    await user.type(input, 'ord-1');
    await user.click(screen.getByRole('button', { name: /add/i }));

    expect(input).toHaveValue('');
  });

  it('adds on Enter', async () => {
    const user = userEvent.setup();
    const { onUpdateItemList, input } = renderList();

    await user.type(input, 'ord-1{Enter}');

    expect(onUpdateItemList).toHaveBeenCalledWith(['ord-1']);
  });

  it('ignores Enter on whitespace', async () => {
    const user = userEvent.setup();
    const { onUpdateItemList, input } = renderList();

    await user.type(input, '   {Enter}');

    expect(onUpdateItemList).not.toHaveBeenCalled();
  });

  it('disables the add button until there is something to add', async () => {
    const user = userEvent.setup();
    const { input } = renderList();
    const addButton = screen.getByRole('button', { name: /add/i });

    expect(addButton).toBeDisabled();
    await user.type(input, '   ');
    expect(addButton).toBeDisabled();
    await user.type(input, 'ord-1');
    expect(addButton).toBeEnabled();
  });

  it('deletes the item at the clicked position', async () => {
    const user = userEvent.setup();
    const { onUpdateItemList } = renderList(['ord-1', 'ord-2', 'ord-3']);

    // The first button in the list is the delete for the first item.
    const deleteButtons = screen.getAllByRole('button');
    await user.click(deleteButtons[1]);

    expect(onUpdateItemList).toHaveBeenCalledWith(['ord-1', 'ord-3']);
    expect(screen.queryByText('ord-2')).not.toBeInTheDocument();
  });

  it('follows a change of the incoming list', () => {
    const onUpdateItemList = vi.fn();
    const { rerender } = render(
      <SearchItemList
        title="Reaction IDs"
        itemList={['ord-1']}
        onUpdateItemList={onUpdateItemList}
      />,
    );
    expect(screen.getByText('ord-1')).toBeInTheDocument();

    rerender(
      <SearchItemList
        title="Reaction IDs"
        itemList={['ord-9']}
        onUpdateItemList={onUpdateItemList}
      />,
    );
    expect(screen.queryByText('ord-1')).not.toBeInTheDocument();
    expect(screen.getByText('ord-9')).toBeInTheDocument();
  });
});
