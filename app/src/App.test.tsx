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

import { act, render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { create, toBinary } from '@bufbuild/protobuf';
import { ReactionSchema } from '@buf/open-reaction-database_ord-schema.bufbuild_es/ord-schema/proto/reaction_pb';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import App from './App';

// Every route is reachable from a plain URL, so the router is exercised by
// pushing the path onto jsdom's history before rendering.
const renderAt = (path: string) => {
  window.history.pushState({}, '', path);
  return render(<App />);
};

// Moves the mounted router to another path, as the browser's back and forward do.
const navigateTo = (path: string) =>
  act(() => {
    window.history.pushState({}, '', path);
    window.dispatchEvent(new PopStateEvent('popstate'));
  });

// A serialized reaction with one empty input per key, base64-encoded the way the
// API returns it.
const encodedReaction = (...inputKeys: string[]): string => {
  const reaction = create(ReactionSchema, {
    inputs: Object.fromEntries(inputKeys.map(key => [key, {}])),
  });
  return btoa(String.fromCharCode(...toBinary(ReactionSchema, reaction)));
};

// Answers /api/reactions from `protos`, keyed by reaction ID. A promise holds the
// response until it resolves.
const stubReactions = (protos: Record<string, string | Promise<string>>) =>
  vi.stubGlobal(
    'fetch',
    vi.fn(async (url: string, init?: RequestInit) => {
      if (url === '/api/reactions') {
        const [reactionId] = JSON.parse(init!.body as string).reaction_ids;
        const proto = await protos[reactionId];
        return { ok: true, status: 200, json: async () => [{ proto }] };
      }
      return { ok: true, status: 200, text: async () => '' };
    }),
  );

const inputTabs = (container: HTMLElement): string[] =>
  [...container.querySelectorAll('#inputs .tab')].map(tab => tab.textContent ?? '');

beforeEach(() => {
  vi.stubGlobal(
    'fetch',
    vi.fn(() => new Promise<Response>(() => {})),
  );
});

afterEach(() => {
  vi.unstubAllGlobals();
  window.history.pushState({}, '', '/');
});

describe('App', () => {
  it('frames every page with the nav and the footer', () => {
    renderAt('/');
    expect(screen.getByRole('link', { name: 'Browse' })).toBeInTheDocument();
    expect(screen.getByRole('link', { name: 'GitHub' })).toBeInTheDocument();
  });

  it('routes / to the home page', () => {
    const { container } = renderAt('/');
    expect(container.querySelector('.home')).toBeInTheDocument();
  });

  it('routes /about to the about page', () => {
    renderAt('/about');
    expect(
      screen.getByRole('heading', { name: 'About', level: 3 }),
    ).toBeInTheDocument();
  });

  it('routes /search to the search page', () => {
    renderAt('/search');
    expect(screen.getByText(/Enter search criteria/)).toBeInTheDocument();
  });

  it('routes /ask to the natural-language search page', () => {
    renderAt('/ask');
    expect(
      screen.getByRole('heading', { name: 'Ask about reactions' }),
    ).toBeInTheDocument();
  });

  it('routes /browse to the dataset list', () => {
    const { container } = renderAt('/browse');
    expect(container.querySelector('#browse-main')).toBeInTheDocument();
  });

  it('routes /selected-set to the reaction set', () => {
    const { container } = renderAt('/selected-set?reaction_id=ord-1');
    expect(container.querySelector('#selected-set-main')).toBeInTheDocument();
  });

  it('routes /dataset/:datasetId to the dataset view', () => {
    renderAt('/dataset/ord_dataset-1');
    expect(screen.getByRole('heading', { name: 'Dataset View' })).toBeInTheDocument();
  });

  it('routes /id/:reactionId to the reaction view', () => {
    const { container } = renderAt('/id/ord-1');
    expect(container.querySelector('.main-reaction-view')).toBeInTheDocument();
  });

  // The selected tabs belong to the reaction they were picked on; the second
  // input tab does not exist on a reaction with one input.
  it('starts a newly routed reaction on its first input', async () => {
    const user = userEvent.setup();
    stubReactions({
      'ord-1': encodedReaction('m1', 'm2'),
      'ord-2': encodedReaction('m3'),
    });
    const { container } = renderAt('/id/ord-1');
    await user.click(await screen.findByText('m2'));

    navigateTo('/id/ord-2');

    expect(await screen.findByText('m3')).toBeInTheDocument();
    expect(inputTabs(container)).toEqual(['m3']);
    expect(container.querySelector('#inputs .tab.selected')?.textContent).toBe('m3');
  });

  it('drops the previous reaction while the next one loads', async () => {
    stubReactions({
      'ord-1': encodedReaction('m1'),
      'ord-2': new Promise(() => {}),
    });
    const { container } = renderAt('/id/ord-1');
    await screen.findByText('m1');

    navigateTo('/id/ord-2');

    expect(screen.queryByText('m1')).toBeNull();
    expect(container.querySelector('.spinner-main')).toBeInTheDocument();
  });

  it('ignores a response for the reaction it routed away from', async () => {
    let respondToFirst!: (proto: string) => void;
    stubReactions({
      'ord-1': new Promise(resolve => (respondToFirst = resolve)),
      'ord-2': encodedReaction('m3'),
    });
    const { container } = renderAt('/id/ord-1');
    navigateTo('/id/ord-2');
    await screen.findByText('m3');

    // A macrotask runs only after the response's whole promise chain has settled.
    await act(async () => {
      respondToFirst(encodedReaction('m1', 'm2'));
      await new Promise(resolve => setTimeout(resolve, 0));
    });

    expect(inputTabs(container)).toEqual(['m3']);
  });

  it('renders no page content for an unknown route', () => {
    const { container } = renderAt('/nope');
    expect(screen.getByRole('link', { name: 'Browse' })).toBeInTheDocument();
    expect(container.querySelector('.home')).toBeNull();
  });
});
