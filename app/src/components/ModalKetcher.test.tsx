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

import { act, fireEvent, render, screen } from '@testing-library/react';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import ModalKetcher from './ModalKetcher';

const setMolecule = vi.fn();
const getSmiles = vi.fn(async () => 'CC=O');

// Ketcher itself never boots under jsdom, so stand in for the API the iframe
// would eventually expose on its contentWindow.
const exposeKetcher = (ketcher: unknown = { setMolecule, getSmiles }) => {
  vi.spyOn(HTMLIFrameElement.prototype, 'contentWindow', 'get').mockReturnValue({
    ketcher,
  } as unknown as Window);
};

const renderModal = (smiles = 'CCO') => {
  const onUpdateSmiles = vi.fn();
  const onCloseModal = vi.fn();
  const result = render(
    <ModalKetcher
      smiles={smiles}
      onUpdateSmiles={onUpdateSmiles}
      onCloseModal={onCloseModal}
    />,
  );
  return { ...result, onUpdateSmiles, onCloseModal };
};

// Drives the fake clock and flushes the state updates and promise
// continuations the elapsed timers kicked off.
const advance = async (ms: number) => {
  await act(async () => {
    await vi.advanceTimersByTimeAsync(ms);
  });
};

// The component polls the iframe once a second until Ketcher shows up.
const waitForKetcher = async () => {
  await advance(1000);
  expect(screen.getByText('Save')).toBeEnabled();
};

beforeEach(() => {
  vi.useFakeTimers();
  setMolecule.mockClear();
  getSmiles.mockClear().mockResolvedValue('CC=O');
  vi.stubGlobal(
    'fetch',
    vi.fn(async () => ({ ok: true, status: 200, json: async () => 'MOLFILE' })),
  );
});

afterEach(() => {
  vi.useRealTimers();
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});

describe('ModalKetcher', () => {
  it('embeds the standalone Ketcher bundle', () => {
    renderModal();
    expect(screen.getByTitle('Ketcher Molecular Editor')).toHaveAttribute(
      'src',
      '/ketcher/index.html',
    );
  });

  it('spins and blocks saving until Ketcher loads', () => {
    const { container } = renderModal();
    expect(container.querySelector('.modal-loading')).toBeInTheDocument();
    expect(screen.getByText('Save')).toBeDisabled();
  });

  it('stops spinning once Ketcher appears', async () => {
    exposeKetcher();
    const { container } = renderModal();

    await waitForKetcher();

    expect(container.querySelector('.modal-loading')).toBeNull();
  });

  // Nothing is going to arrive after 30s; give up rather than spin forever.
  it('gives up after thirty seconds', async () => {
    const consoleWarn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    const { container } = renderModal();

    await advance(30_000);

    expect(container.querySelector('.modal-loading')).toBeNull();
    expect(consoleWarn).toHaveBeenCalledWith(
      'Ketcher failed to load within 30 seconds',
    );
  });

  describe('loading the starting structure', () => {
    it('converts the SMILES to a molfile and hands it to Ketcher', async () => {
      exposeKetcher();
      renderModal('CCO');

      await waitForKetcher();

      expect(fetch).toHaveBeenCalledWith('/api/molfile?smiles=CCO');
      expect(setMolecule).toHaveBeenCalledWith('MOLFILE');
    });

    it('url-encodes the SMILES', async () => {
      exposeKetcher();
      renderModal('C/C=C/C');

      await waitForKetcher();

      expect(fetch).toHaveBeenCalledWith('/api/molfile?smiles=C%2FC%3DC%2FC');
    });

    it('opens an empty editor when there is no starting structure', async () => {
      exposeKetcher();
      renderModal('');

      await waitForKetcher();

      expect(fetch).not.toHaveBeenCalled();
      expect(setMolecule).not.toHaveBeenCalled();
    });

    it('warns rather than drawing when the conversion fails', async () => {
      const consoleWarn = vi.spyOn(console, 'warn').mockImplementation(() => {});
      exposeKetcher();
      vi.stubGlobal(
        'fetch',
        vi.fn(async () => ({ ok: false, status: 500 })),
      );
      renderModal('CCO');

      await waitForKetcher();

      expect(consoleWarn).toHaveBeenCalled();
      expect(setMolecule).not.toHaveBeenCalled();
    });

    it('survives a network error', async () => {
      const consoleError = vi.spyOn(console, 'error').mockImplementation(() => {});
      exposeKetcher();
      vi.stubGlobal(
        'fetch',
        vi.fn(async () => Promise.reject(new Error('offline'))),
      );
      renderModal('CCO');

      await waitForKetcher();

      expect(consoleError).toHaveBeenCalled();
    });
  });

  describe('saving', () => {
    it('reports the drawn structure and closes', async () => {
      exposeKetcher();
      const { onUpdateSmiles, onCloseModal } = renderModal();
      await waitForKetcher();

      await act(async () => {
        fireEvent.click(screen.getByText('Save'));
      });

      expect(onUpdateSmiles).toHaveBeenCalledWith('CC=O');
      expect(onCloseModal).toHaveBeenCalled();
    });

    it('stays open when Ketcher cannot produce a SMILES', async () => {
      const consoleError = vi.spyOn(console, 'error').mockImplementation(() => {});
      getSmiles.mockRejectedValue(new Error('no structure'));
      exposeKetcher();
      const { onUpdateSmiles, onCloseModal } = renderModal();
      await waitForKetcher();

      await act(async () => {
        fireEvent.click(screen.getByText('Save'));
      });

      expect(consoleError).toHaveBeenCalled();
      expect(onUpdateSmiles).not.toHaveBeenCalled();
      expect(onCloseModal).not.toHaveBeenCalled();
    });
  });

  describe('dismissing', () => {
    it('closes on Cancel', () => {
      const { onCloseModal } = renderModal();
      fireEvent.click(screen.getByText('Cancel'));
      expect(onCloseModal).toHaveBeenCalledTimes(1);
    });

    it('closes on a click outside the editor', () => {
      const { container, onCloseModal } = renderModal();
      fireEvent.click(container.querySelector('.background')!);
      expect(onCloseModal).toHaveBeenCalledTimes(1);
    });

    it('stays open on a click inside the editor', () => {
      const { container, onCloseModal } = renderModal();
      fireEvent.click(container.querySelector('#ketcher_modal')!);
      expect(onCloseModal).not.toHaveBeenCalled();
    });

    it('closes on Escape', () => {
      const { onCloseModal } = renderModal();
      fireEvent.keyDown(document, { key: 'Escape' });
      expect(onCloseModal).toHaveBeenCalledTimes(1);
    });

    it('ignores other keys', () => {
      const { onCloseModal } = renderModal();
      fireEvent.keyDown(document, { key: 'Enter' });
      expect(onCloseModal).not.toHaveBeenCalled();
    });

    // The listeners and the poll must not outlive the modal.
    it('stops polling once unmounted', async () => {
      exposeKetcher();
      const { unmount, onCloseModal } = renderModal();

      unmount();
      await advance(30_000);
      fireEvent.keyDown(document, { key: 'Escape' });

      expect(onCloseModal).not.toHaveBeenCalled();
    });
  });
});
