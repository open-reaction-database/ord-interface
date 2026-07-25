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

import { render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { afterEach, describe, expect, it, vi } from 'vitest';
import CopyButton from './CopyButton';

// userEvent.setup() installs its own navigator.clipboard stub, so spy on that
// rather than replacing the clipboard ourselves.
const setup = () => {
  const user = userEvent.setup();
  return { user, writeText: vi.spyOn(navigator.clipboard, 'writeText') };
};

afterEach(() => {
  vi.restoreAllMocks();
  vi.useRealTimers();
});

describe('CopyButton', () => {
  it('copies the text to the clipboard', async () => {
    const { user, writeText } = setup();
    render(<CopyButton textToCopy="CC(=O)O" />);

    await user.click(screen.getByRole('button'));

    expect(writeText).toHaveBeenCalledWith('CC(=O)O');
    await expect(navigator.clipboard.readText()).resolves.toBe('CC(=O)O');
  });

  it('confirms the copy, then withdraws the confirmation', async () => {
    const { user } = setup();
    render(<CopyButton textToCopy="CC(=O)O" />);

    await user.click(screen.getByRole('button'));
    expect(await screen.findByText('Copied to clipboard!')).toBeInTheDocument();

    await waitFor(
      () => expect(screen.queryByText('Copied to clipboard!')).not.toBeInTheDocument(),
      { timeout: 3000 },
    );
  });

  it('shows the default icon and no caption', () => {
    render(<CopyButton textToCopy="CC(=O)O" />);
    expect(screen.getByText('content_copy')).toBeInTheDocument();
    expect(screen.getByRole('button').querySelector('.copy')).toBeNull();
  });

  it('shows a custom icon and caption', () => {
    render(
      <CopyButton
        textToCopy="CC(=O)O"
        icon="share"
        buttonText="Shareable Link"
      />,
    );
    expect(screen.getByText('share')).toBeInTheDocument();
    expect(screen.getByText('Shareable Link')).toBeInTheDocument();
  });

  // A denied clipboard permission must not claim the copy succeeded.
  it('does not confirm when the clipboard write fails', async () => {
    const { user, writeText } = setup();
    writeText.mockRejectedValue(new Error('denied'));
    const consoleError = vi.spyOn(console, 'error').mockImplementation(() => {});
    render(<CopyButton textToCopy="CC(=O)O" />);

    await user.click(screen.getByRole('button'));

    await waitFor(() => expect(consoleError).toHaveBeenCalled());
    expect(screen.queryByText('Copied to clipboard!')).not.toBeInTheDocument();
  });
});
