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

import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { renderHook, waitFor } from '@testing-library/react';
import { createElement, type ReactNode } from 'react';
import reaction_pb from 'ord-schema';
import { afterEach, describe, expect, it, vi } from 'vitest';
import type { NLQueryResponse } from '../types/search';
import { useNLQuery } from './useNLQuery';

const encodedReaction = (reactionId: string): string => {
  const reaction = new reaction_pb.Reaction();
  reaction.setReactionId(reactionId);
  return btoa(String.fromCharCode(...reaction.serializeBinary()));
};

const wrapper = ({ children }: { children: ReactNode }) =>
  createElement(
    QueryClientProvider,
    {
      client: new QueryClient({
        defaultOptions: { queries: { retry: false, gcTime: 0 } },
      }),
    },
    children,
  );

const response = (overrides: Partial<NLQueryResponse> = {}): NLQueryResponse => ({
  query: 'suzuki couplings with high yield',
  interpretation: { components: [], min_yield: 80 },
  resolved_components: [
    {
      identifier: 'phenylboronic acid',
      smiles: 'OB(O)c1ccccc1',
      resolver: 'pubchem',
      target: 'INPUT',
      mode: 'EXACT',
    },
  ],
  query_components: ['{"pattern":"OB(O)c1ccccc1","target":"input","mode":"exact"}'],
  results: [{ reaction_id: 'ord-1', proto: encodedReaction('ord-1') }],
  dry_run: false,
  ...overrides,
});

const stubFetch = (body: NLQueryResponse, status = 200) => {
  const fetchMock = vi.fn(
    async () =>
      ({ ok: status < 400, status, json: async () => body }) as unknown as Response,
  );
  vi.stubGlobal('fetch', fetchMock);
  return fetchMock;
};

const renderNLQuery = (query: string | null, enabled = true, dryRun = false) =>
  renderHook(() => useNLQuery(query, enabled, dryRun), { wrapper });

afterEach(() => {
  vi.unstubAllGlobals();
});

describe('useNLQuery', () => {
  it('stays idle when disabled', () => {
    const fetchMock = stubFetch(response());
    renderNLQuery('suzuki couplings', false);
    expect(fetchMock).not.toHaveBeenCalled();
  });

  it('stays idle without a question', () => {
    const fetchMock = stubFetch(response());
    renderNLQuery(null);
    expect(fetchMock).not.toHaveBeenCalled();
  });

  // A box the user has only typed spaces into is not a question.
  it('stays idle for a blank question', () => {
    const fetchMock = stubFetch(response());
    renderNLQuery('   ');
    expect(fetchMock).not.toHaveBeenCalled();
  });

  it('url-encodes the question', async () => {
    const fetchMock = stubFetch(response());
    const { result } = renderNLQuery('yield > 80% & Pd');

    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(fetchMock).toHaveBeenCalledWith(
      '/api/nl_query?q=yield%20%3E%2080%25%20%26%20Pd',
      undefined,
    );
  });

  it('asks for a dry run when requested', async () => {
    const fetchMock = stubFetch(response({ dry_run: true }));
    const { result } = renderNLQuery('suzuki couplings', true, true);

    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(fetchMock).toHaveBeenCalledWith(
      '/api/nl_query?q=suzuki%20couplings&dry_run=true',
      undefined,
    );
    expect(result.current.data?.dryRun).toBe(true);
  });

  it('renames the payload fields for the view layer', async () => {
    const payload = response();
    stubFetch(payload);
    const { result } = renderNLQuery('suzuki couplings');

    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(result.current.data?.interpretation).toEqual(payload.interpretation);
    expect(result.current.data?.resolvedComponents).toEqual(
      payload.resolved_components,
    );
    expect(result.current.data?.queryComponents).toEqual(payload.query_components);
    expect(result.current.data?.dryRun).toBe(false);
  });

  it('deserializes the result protos', async () => {
    stubFetch(response());
    const { result } = renderNLQuery('suzuki couplings');

    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(result.current.data?.results).toEqual([
      expect.objectContaining({
        reaction_id: 'ord-1',
        data: expect.objectContaining({ reactionId: 'ord-1' }),
      }),
    ]);
  });

  it('surfaces a failed request', async () => {
    stubFetch(response(), 500);
    const { result } = renderNLQuery('suzuki couplings');

    await waitFor(() => expect(result.current.isError).toBe(true));
    expect(result.current.error?.message).toBe('nl_query failed (HTTP 500)');
  });

  // Dry runs and real runs are cached separately, so toggling the switch
  // re-queries instead of replaying the other mode's answer.
  it('re-queries when the dry-run flag changes', async () => {
    const fetchMock = stubFetch(response());
    const { result, rerender } = renderHook(
      ({ dryRun }: { dryRun: boolean }) => useNLQuery('suzuki couplings', true, dryRun),
      { wrapper, initialProps: { dryRun: false } },
    );
    await waitFor(() => expect(result.current.isSuccess).toBe(true));

    rerender({ dryRun: true });
    await waitFor(() => expect(fetchMock).toHaveBeenCalledTimes(2));
  });
});
