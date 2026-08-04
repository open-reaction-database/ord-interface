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
import { useSearchTask } from './useSearchTask';

// A serialized Reaction, base64-encoded the way the API returns it.
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

const renderSearchTask = (queryString: string | null, enabled = true) =>
  renderHook(() => useSearchTask(queryString, enabled), { wrapper });

const jsonResponse = (body: unknown, status = 200) =>
  ({ ok: status < 400, status, json: async () => body }) as Response;

// Stubs the two-step protocol: submit_query hands back a task id, and
// fetch_query_result replays the given statuses in order.
const stubProtocol = (
  results: Array<{ status: number; body?: unknown }>,
  taskId = 'task-1',
) => {
  let poll = 0;
  const fetchMock = vi.fn(async (url: string) => {
    if (url.startsWith('/api/submit_query')) return jsonResponse(taskId);
    const next = results[Math.min(poll, results.length - 1)];
    poll += 1;
    return jsonResponse(next.body ?? [], next.status);
  });
  vi.stubGlobal('fetch', fetchMock);
  return fetchMock;
};

const submitCalls = (fetchMock: ReturnType<typeof vi.fn>): string[] =>
  fetchMock.mock.calls
    .map(call => call[0] as string)
    .filter(url => url.startsWith('/api/submit_query'));

afterEach(() => {
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});

describe('useSearchTask', () => {
  it('stays idle when disabled', () => {
    const fetchMock = stubProtocol([{ status: 200 }]);
    renderSearchTask('?dataset_id=ord_dataset-1', false);
    expect(fetchMock).not.toHaveBeenCalled();
  });

  it('stays idle without a query string', () => {
    const fetchMock = stubProtocol([{ status: 200 }]);
    renderSearchTask(null);
    expect(fetchMock).not.toHaveBeenCalled();
  });

  it('submits the query string verbatim', async () => {
    const fetchMock = stubProtocol([{ status: 200 }]);
    const { result } = renderSearchTask('?dataset_id=ord_dataset-1');

    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(submitCalls(fetchMock)).toEqual([
      '/api/submit_query?dataset_id=ord_dataset-1',
    ]);
  });

  it('polls the task id returned by submit_query', async () => {
    const fetchMock = stubProtocol([{ status: 200 }], 'task-42');
    const { result } = renderSearchTask('?dataset_id=ord_dataset-1');

    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(fetchMock).toHaveBeenCalledWith('/api/fetch_query_result?task_id=task-42');
  });

  it('deserializes the result protos', async () => {
    stubProtocol([
      {
        status: 200,
        body: [
          {
            reaction_id: 'ord-1',
            dataset_id: 'ord_dataset-1',
            proto: encodedReaction('ord-1'),
          },
        ],
      },
    ]);
    const { result } = renderSearchTask('?dataset_id=ord_dataset-1');

    await waitFor(() => expect(result.current.isSuccess).toBe(true));
    expect(result.current.data).toEqual({
      status: 'success',
      results: [
        expect.objectContaining({
          reaction_id: 'ord-1',
          dataset_id: 'ord_dataset-1',
          data: expect.objectContaining({ reactionId: 'ord-1' }),
        }),
      ],
    });
  });

  it('reports pending while the task is still running', async () => {
    stubProtocol([{ status: 202 }]);
    const { result } = renderSearchTask('?dataset_id=ord_dataset-1');

    await waitFor(() =>
      expect(result.current.data).toEqual({ status: 'pending', taskId: 'task-1' }),
    );
  });

  it('keeps polling until the task completes', async () => {
    stubProtocol([
      { status: 202 },
      { status: 202 },
      {
        status: 200,
        body: [{ reaction_id: 'ord-1', proto: encodedReaction('ord-1') }],
      },
    ]);
    const { result } = renderSearchTask('?dataset_id=ord_dataset-1');

    await waitFor(() => expect(result.current.data?.status).toBe('success'), {
      timeout: 5000,
    });
  });

  // Submitting twice would leave an orphaned task running on the backend.
  it('submits only once across the whole poll cycle', async () => {
    const fetchMock = stubProtocol([{ status: 202 }, { status: 202 }, { status: 200 }]);
    const { result } = renderSearchTask('?dataset_id=ord_dataset-1');

    await waitFor(() => expect(result.current.data?.status).toBe('success'), {
      timeout: 5000,
    });
    expect(submitCalls(fetchMock)).toHaveLength(1);
  });

  it('surfaces a failed submit', async () => {
    vi.stubGlobal(
      'fetch',
      vi.fn(async () => jsonResponse({}, 500)),
    );
    const { result } = renderSearchTask('?dataset_id=ord_dataset-1');

    await waitFor(() => expect(result.current.isError).toBe(true));
    expect(result.current.error?.message).toBe('submit_query failed (HTTP 500)');
  });

  it('surfaces a failed poll', async () => {
    stubProtocol([{ status: 500 }], 'task-7');
    const { result } = renderSearchTask('?dataset_id=ord_dataset-1');

    await waitFor(() => expect(result.current.isError).toBe(true));
    expect(result.current.error?.message).toBe('Search task task-7 failed (HTTP 500)');
  });

  it('gives up on a task that runs past the timeout', async () => {
    let now = 0;
    vi.spyOn(Date, 'now').mockImplementation(() => now);
    vi.stubGlobal(
      'fetch',
      vi.fn(async (url: string) => {
        if (url.startsWith('/api/submit_query')) {
          // The clock passes the deadline while the submit is in flight.
          now = 200_000;
          return jsonResponse('task-9');
        }
        return jsonResponse([], 200);
      }),
    );
    const { result } = renderSearchTask('?dataset_id=ord_dataset-1');

    await waitFor(() => expect(result.current.isError).toBe(true));
    expect(result.current.error?.message).toBe(
      'Search task task-9 timed out after 120s',
    );
  });

  it('starts a fresh task when the query string changes', async () => {
    const fetchMock = stubProtocol([{ status: 200 }]);
    const { result, rerender } = renderHook(
      ({ query }: { query: string }) => useSearchTask(query, true),
      { wrapper, initialProps: { query: '?dataset_id=ord_dataset-1' } },
    );
    await waitFor(() => expect(result.current.isSuccess).toBe(true));

    rerender({ query: '?dataset_id=ord_dataset-2' });
    await waitFor(() => expect(result.current.isSuccess).toBe(true));

    expect(submitCalls(fetchMock)).toEqual([
      '/api/submit_query?dataset_id=ord_dataset-1',
      '/api/submit_query?dataset_id=ord_dataset-2',
    ]);
  });
});
