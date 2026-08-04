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
import reaction_pb from 'ord-schema';
import { MemoryRouter } from 'react-router-dom';
import { afterEach, describe, expect, it, vi } from 'vitest';
import MainSelectedSet from './MainSelectedSet';

const encodedReaction = (reactionId: string): string => {
  const reaction = new reaction_pb.Reaction();
  reaction.setReactionId(reactionId);
  return btoa(String.fromCharCode(...reaction.serializeBinary()));
};

// The bodies POSTed to /api/reactions, for the assertions below.
const posted: RequestInit[] = [];

// /api/reactions returns the set; /api/reaction_summary draws each card.
const stubApi = (reactionIds: string[], status = 200) => {
  const fetchMock = vi.fn(async (url: string, init?: RequestInit) => {
    if (url === '/api/reactions') {
      posted.push(init!);
      return {
        ok: status < 400,
        status,
        json: async () =>
          reactionIds.map(id => ({
            reaction_id: id,
            dataset_id: 'ord_dataset-1',
            proto: encodedReaction(id),
          })),
      };
    }
    return { ok: true, status: 200, text: async () => '' };
  });
  vi.stubGlobal('fetch', fetchMock);
  return fetchMock;
};

const renderSet = (search: string) =>
  render(
    <MemoryRouter initialEntries={[`/selected-set${search}`]}>
      <MainSelectedSet />
    </MemoryRouter>,
  );

afterEach(() => {
  posted.length = 0;
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});

describe('MainSelectedSet', () => {
  it('posts the reaction ids from the URL', async () => {
    stubApi(['ord-1', 'ord-2']);
    renderSet('?reaction_id=ord-1&reaction_id=ord-2');

    expect(await screen.findByText('Reaction Set')).toBeInTheDocument();
    expect(JSON.parse(posted[0].body as string)).toEqual({
      reaction_ids: ['ord-1', 'ord-2'],
    });
  });

  it('renders a card per reaction', async () => {
    stubApi(['ord-1', 'ord-2']);
    const { container } = renderSet('?reaction_id=ord-1&reaction_id=ord-2');

    await screen.findByText('Reaction Set');
    expect(container.querySelectorAll('.reaction-container')).toHaveLength(2);
  });

  it('offers a shareable link and a download', async () => {
    stubApi(['ord-1']);
    renderSet('?reaction_id=ord-1');

    expect(await screen.findByText('Shareable Link')).toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Download Reaction Set' })).toBeEnabled();
  });

  it('opens the download modal', async () => {
    const user = userEvent.setup();
    stubApi(['ord-1']);
    renderSet('?reaction_id=ord-1');

    await user.click(
      await screen.findByRole('button', { name: 'Download Reaction Set' }),
    );

    expect(screen.getByText('Download Results')).toBeInTheDocument();
  });

  // Landing here with no ids in the URL should settle on the empty state
  // rather than spinning forever.
  it('stops loading when the URL carries no ids', async () => {
    const fetchMock = stubApi([]);
    renderSet('');

    expect(
      await screen.findByText('There was an issue fetching your selected reactions.'),
    ).toBeInTheDocument();
    expect(fetchMock).not.toHaveBeenCalled();
  });

  it('reports a failed fetch instead of an empty page', async () => {
    const consoleError = vi.spyOn(console, 'error').mockImplementation(() => {});
    stubApi(['ord-1'], 500);
    renderSet('?reaction_id=ord-1');

    expect(
      await screen.findByText('There was an issue fetching your selected reactions.'),
    ).toBeInTheDocument();
    expect(consoleError).toHaveBeenCalled();
  });

  it('spins while the reactions are in flight', () => {
    vi.stubGlobal(
      'fetch',
      vi.fn(() => new Promise<Response>(() => {})),
    );
    const { container } = renderSet('?reaction_id=ord-1');
    expect(container.querySelector('.spinner-main')).toBeInTheDocument();
  });
});
