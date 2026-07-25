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
import FloatingModal from './FloatingModal';

const renderModal = (className?: string) => {
  const onCloseModal = vi.fn();
  const { container } = render(
    <FloatingModal
      title="Download Results"
      onCloseModal={onCloseModal}
      className={className}
    >
      <p>body content</p>
    </FloatingModal>,
  );
  return { container, onCloseModal };
};

describe('FloatingModal', () => {
  it('shows the title and the children', () => {
    renderModal();
    expect(screen.getByText('Download Results')).toBeInTheDocument();
    expect(screen.getByText('body content')).toBeInTheDocument();
  });

  it('closes on the close control', async () => {
    const user = userEvent.setup();
    const { container, onCloseModal } = renderModal();

    await user.click(container.querySelector('.close')!);

    expect(onCloseModal).toHaveBeenCalledTimes(1);
  });

  it('closes on a click outside the modal', async () => {
    const user = userEvent.setup();
    const { container, onCloseModal } = renderModal();

    await user.click(container.querySelector('.modal-background')!);

    expect(onCloseModal).toHaveBeenCalledTimes(1);
  });

  it('stays open on a click inside the modal', async () => {
    const user = userEvent.setup();
    const { onCloseModal } = renderModal();

    await user.click(screen.getByText('body content'));

    expect(onCloseModal).not.toHaveBeenCalled();
  });

  it('appends a caller-supplied class without dropping its own', () => {
    const { container } = renderModal('wide');
    expect(container.querySelector('.modal-container')).toHaveClass(
      'modal-container',
      'wide',
    );
  });

  it('leaves the class list clean when no extra class is given', () => {
    const { container } = renderModal();
    expect(container.querySelector('.modal-container')).toHaveAttribute(
      'class',
      'modal-container',
    );
  });
});
