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
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import CompoundView from './CompoundView';

const SMILES = 2;
const NAME = 6;
const REAGENT_ROLE = 2;

const renderCompound = (component: unknown) =>
  render(<CompoundView component={component as never} />);

const openRawData = async (user: ReturnType<typeof userEvent.setup>) => {
  await user.click(screen.getByText('<>'));
};

// Each raw-data line is a <pre>, and JSX splits the label and the JSON into
// separate text nodes, so read the lines whole rather than by text match.
const rawLines = (container: HTMLElement): string[] =>
  [...container.querySelectorAll('.data pre')].map(line => line.textContent ?? '');

beforeEach(() => {
  vi.stubGlobal(
    'fetch',
    vi.fn(async () => ({
      ok: true,
      status: 200,
      json: async () => '<svg class="mol"></svg>',
    })),
  );
});

afterEach(() => {
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});

describe('CompoundView', () => {
  it('renders without a component', () => {
    const { container } = renderCompound(undefined);
    expect(container.querySelector('.compound-view')).toBeInTheDocument();
    expect(fetch).not.toHaveBeenCalled();
  });

  it('shows the amount column only when there is an amount', () => {
    const { unmount } = renderCompound({});
    expect(screen.queryByText('Amount')).not.toBeInTheDocument();
    unmount();

    renderCompound({ amount: { mass: { value: 250, units: 3 } } });
    expect(screen.getByText('Amount')).toBeInTheDocument();
    expect(screen.getByText('250 milligram')).toBeInTheDocument();
  });

  it('shows the reaction role in lowercase', () => {
    const { container } = renderCompound({ reactionRole: REAGENT_ROLE });
    expect(container.querySelector('.role')).toHaveTextContent('reagent');
  });

  it('leaves the role blank when none was recorded', () => {
    const { container } = renderCompound({});
    expect(container.querySelector('.role')).toBeEmptyDOMElement();
  });

  describe('the compound drawing', () => {
    it('posts the SMILES identifier and renders the returned SVG', async () => {
      const { container } = renderCompound({
        identifiersList: [{ type: SMILES, value: 'CCO' }],
      });

      await waitFor(() =>
        expect(container.querySelector('.svg')?.innerHTML).toBe(
          '<svg class="mol"></svg>',
        ),
      );
      expect(fetch).toHaveBeenCalledWith(
        '/api/compound_svg',
        expect.objectContaining({
          method: 'POST',
          headers: { 'Content-Type': 'application/x-protobuf' },
        }),
      );
    });

    // Without a SMILES there is nothing to draw, so no request should go out.
    it('skips the request when no identifier is a SMILES', () => {
      renderCompound({ identifiersList: [{ type: NAME, value: 'ethanol' }] });
      expect(fetch).not.toHaveBeenCalled();
    });

    // The 4xx/5xx body is an HTML error page, not an SVG.
    it('does not render an error page in the drawing slot', async () => {
      const consoleError = vi.spyOn(console, 'error').mockImplementation(() => {});
      const json = vi.fn();
      vi.stubGlobal(
        'fetch',
        vi.fn(async () => ({ ok: false, status: 500, json })),
      );

      const { container } = renderCompound({
        identifiersList: [{ type: SMILES, value: 'CCO' }],
      });

      await waitFor(() => expect(consoleError).toHaveBeenCalled());
      expect(json).not.toHaveBeenCalled();
      expect(container.querySelector('.svg')).toBeEmptyDOMElement();
    });

    it('survives a network error', async () => {
      const consoleError = vi.spyOn(console, 'error').mockImplementation(() => {});
      vi.stubGlobal(
        'fetch',
        vi.fn(async () => Promise.reject(new Error('offline'))),
      );

      const { container } = renderCompound({
        identifiersList: [{ type: SMILES, value: 'CCO' }],
      });

      await waitFor(() => expect(consoleError).toHaveBeenCalled());
      expect(container.querySelector('.svg')).toBeEmptyDOMElement();
    });
  });

  describe('the raw data modal', () => {
    it('names the identifier types', async () => {
      const user = userEvent.setup();
      const { container } = renderCompound({
        identifiersList: [{ type: SMILES, value: 'CCO' }],
      });

      await openRawData(user);

      expect(rawLines(container)).toContain(
        'identifiers: {"type":"SMILES","value":"CCO"}',
      );
    });

    it('reports the amount under its category', async () => {
      const user = userEvent.setup();
      const { container } = renderCompound({
        amount: { volume: { value: 10, units: 2 } },
      });

      await openRawData(user);

      expect(rawLines(container)).toContain(
        'amount: {"volume":{"value":10,"units":"MILLILITER"}}',
      );
    });

    it('omits the amount when none was recorded', async () => {
      const user = userEvent.setup();
      const { container } = renderCompound({ reactionRole: REAGENT_ROLE });

      await openRawData(user);

      expect(rawLines(container)).toEqual(['reaction_role: REAGENT']);
    });

    it('names the preparation types', async () => {
      const user = userEvent.setup();
      const { container } = renderCompound({
        preparationsList: [{ type: 2, details: 'used as received' }],
      });

      await openRawData(user);

      expect(rawLines(container)).toContain(
        'preparations: {"type":"NONE","details":"used as received"}',
      );
    });

    it('reports the product-only fields', async () => {
      const user = userEvent.setup();
      const { container } = renderCompound({
        isolatedColor: 'white',
        texture: { type: 2, details: 'fine' },
      });

      await openRawData(user);

      expect(rawLines(container)).toContain('isolated_color: white');
      expect(rawLines(container)).toContain('texture: {"type":"2","details":"fine"}');
    });

    it('closes again', async () => {
      const user = userEvent.setup();
      renderCompound({ reactionRole: REAGENT_ROLE });

      await openRawData(user);
      await user.click(screen.getByText('✕'));

      expect(screen.queryByText('Raw Data')).not.toBeInTheDocument();
    });
  });
});
