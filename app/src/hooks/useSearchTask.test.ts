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
import { act, renderHook, waitFor } from '@testing-library/react';
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

const requestedUrls = (fetchMock: ReturnType<typeof vi.fn>): string[] =>
  fetchMock.mock.calls.map(call => call[0] as string);

const submitCalls = (fetchMock: ReturnType<typeof vi.fn>): string[] =>
  requestedUrls(fetchMock).filter(url => url.startsWith('/api/submit_query'));

const reactionIds = (data: unknown): string[] | undefined =>
  (data as { results?: Array<{ reaction_id: string }> } | undefined)?.results?.map(
    result => result.reaction_id,
  );

// Renders the hook against one query client for the whole test, so switching back
// to an earlier query string finds whatever the cache kept for it.
const renderWithSharedClient = (initialQuery: string) => {
  const client = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return renderHook(({ query }: { query: string }) => useSearchTask(query, true), {
    wrapper: ({ children }: { children: ReactNode }) =>
      createElement(QueryClientProvider, { client }, children),
    initialProps: { query: initialQuery },
  });
};

// Starts search A, switches to search B while A's submit_query is still in flight,
// and lets A's submit return once B's task is running. Task A completes on its
// first poll; task B reports pending once, then completes.
const raceSearches = async () => {
  let releaseSubmitA!: () => void;
  let taskBPolls = 0;
  const fetchMock = vi.fn(async (url: string) => {
    if (url === '/api/submit_query?q=A') {
      await new Promise<void>(resolve => (releaseSubmitA = resolve));
      return jsonResponse('task-A');
    }
    if (url === '/api/submit_query?q=B') return jsonResponse('task-B');
    if (url === '/api/fetch_query_result?task_id=task-A') {
      return jsonResponse([{ reaction_id: 'A', proto: encodedReaction('A') }]);
    }
    if (url === '/api/fetch_query_result?task_id=task-B') {
      taskBPolls += 1;
      return taskBPolls === 1
        ? jsonResponse([], 202)
        : jsonResponse([{ reaction_id: 'B', proto: encodedReaction('B') }]);
    }
    return jsonResponse({}, 404);
  });
  vi.stubGlobal('fetch', fetchMock);

  const hook = renderWithSharedClient('?q=A');
  hook.rerender({ query: '?q=B' });
  await waitFor(() =>
    expect(hook.result.current.data).toEqual({ status: 'pending', taskId: 'task-B' }),
  );
  // A macrotask runs only after A's whole promise chain has settled.
  await act(async () => {
    releaseSubmitA();
    await new Promise(resolve => setTimeout(resolve, 0));
  });
  return { fetchMock, ...hook };
};

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
    expect(requestedUrls(fetchMock)).toContain(
      '/api/fetch_query_result?task_id=task-42',
    );
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
        return jsonResponse([], 202);
      }),
    );
    const { result } = renderSearchTask('?dataset_id=ord_dataset-1');

    await waitFor(() => expect(result.current.isError).toBe(true));
    expect(result.current.error?.message).toBe(
      'Search task task-9 timed out after 120s',
    );
  });

  // Polling pauses in a hidden tab, so the next poll can come after the deadline.
  it('returns a result that is ready once the deadline has passed', async () => {
    let now = 0;
    vi.spyOn(Date, 'now').mockImplementation(() => now);
    stubProtocol([
      { status: 202 },
      {
        status: 200,
        body: [{ reaction_id: 'ord-1', proto: encodedReaction('ord-1') }],
      },
    ]);
    const { result } = renderSearchTask('?dataset_id=ord_dataset-1');
    await waitFor(() => expect(result.current.data?.status).toBe('pending'));

    now = 200_000;

    await waitFor(() => expect(reactionIds(result.current.data)).toEqual(['ord-1']), {
      timeout: 3000,
    });
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

  // A request with no bound would leave the search loading for as long as it hangs.
  it('fails a search whose poll request stalls', async () => {
    const stall = new AbortController();
    const timeout = vi.spyOn(AbortSignal, 'timeout').mockReturnValue(stall.signal);
    vi.stubGlobal(
      'fetch',
      vi.fn((url: string, init?: RequestInit) =>
        url.startsWith('/api/submit_query')
          ? Promise.resolve(jsonResponse('task-1'))
          : new Promise<Response>((_, reject) =>
              init?.signal?.addEventListener('abort', () =>
                reject(init.signal!.reason),
              ),
            ),
      ),
    );
    const { result } = renderSearchTask('?dataset_id=ord_dataset-1');
    await waitFor(() => expect(timeout).toHaveBeenCalledTimes(2));

    stall.abort(new DOMException('signal timed out', 'TimeoutError'));

    await waitFor(() => expect(result.current.isError).toBe(true));
    expect(timeout).toHaveBeenCalledWith(90_000);
  });

  // Rerunning must not re-read a task whose result could not be decoded.
  it('submits again when a search whose result did not decode runs again', async () => {
    const fetchMock = stubProtocol([{ status: 200, body: [{ proto: '!' }] }]);
    const { result } = renderSearchTask('?dataset_id=ord_dataset-1');
    await waitFor(() => expect(result.current.isError).toBe(true));

    await act(() => result.current.refetch());

    expect(submitCalls(fetchMock)).toHaveLength(2);
  });

  // Asking the overdue task again would time out at once, against its old start.
  it('submits again when a timed-out search runs again', async () => {
    let now = 0;
    vi.spyOn(Date, 'now').mockImplementation(() => now);
    const fetchMock = stubProtocol([{ status: 202 }]);
    const { result } = renderSearchTask('?dataset_id=ord_dataset-1');
    await waitFor(() => expect(result.current.data?.status).toBe('pending'));
    now = 200_000;
    await waitFor(() => expect(result.current.isError).toBe(true), { timeout: 3000 });

    await act(() => result.current.refetch());

    await waitFor(() => expect(result.current.isError).toBe(false));
    expect(submitCalls(fetchMock)).toHaveLength(2);
  });

  // The task may have expired, or its result may be too slow to build again, so a
  // rerun after any failure starts over.
  it.each([
    ['fails in transit', () => Promise.reject(new TypeError('Failed to fetch'))],
    ['gets a server error', () => Promise.resolve(jsonResponse({}, 504))],
  ])('submits again when a search whose poll %s runs again', async (_, failedPoll) => {
    let polls = 0;
    const fetchMock = vi.fn(async (url: string) => {
      if (url.startsWith('/api/submit_query')) return jsonResponse('task-1');
      polls += 1;
      return polls === 1
        ? failedPoll()
        : jsonResponse([{ reaction_id: 'ord-1', proto: encodedReaction('ord-1') }]);
    });
    vi.stubGlobal('fetch', fetchMock);
    const { result } = renderSearchTask('?dataset_id=ord_dataset-1');
    await waitFor(() => expect(result.current.isError).toBe(true));

    await act(() => result.current.refetch());

    await waitFor(() => expect(reactionIds(result.current.data)).toEqual(['ord-1']));
    expect(submitCalls(fetchMock)).toHaveLength(2);
  });

  // The cached result still says pending after the error, and polling on from it
  // would resubmit the query and start another backend task.
  describe('after giving up on a pending task', () => {
    const outlastPollInterval = () => new Promise(resolve => setTimeout(resolve, 1500));

    it('stops polling once the task times out', async () => {
      let now = 0;
      vi.spyOn(Date, 'now').mockImplementation(() => now);
      const fetchMock = stubProtocol([{ status: 202 }]);
      const { result } = renderSearchTask('?dataset_id=ord_dataset-1');
      await waitFor(() => expect(result.current.data?.status).toBe('pending'));

      now = 200_000;
      await waitFor(() => expect(result.current.isError).toBe(true), { timeout: 3000 });
      await outlastPollInterval();

      expect(result.current.isError).toBe(true);
      expect(submitCalls(fetchMock)).toHaveLength(1);
    });

    it('stops polling once a poll fails', async () => {
      const fetchMock = stubProtocol([{ status: 202 }, { status: 500 }]);
      const { result } = renderSearchTask('?dataset_id=ord_dataset-1');
      await waitFor(() => expect(result.current.isError).toBe(true), { timeout: 3000 });
      await outlastPollInterval();

      expect(result.current.isError).toBe(true);
      expect(submitCalls(fetchMock)).toHaveLength(1);
    });
  });

  // React Query marks a query invalidated when a fetch fails, so a failed search
  // runs again on return even though staleTime is Infinity.
  it('runs a failed search again when it is revisited', async () => {
    let submits = 0;
    const fetchMock = vi.fn(async (url: string) => {
      if (url === '/api/submit_query?q=A') return jsonResponse(`task-A${++submits}`);
      if (url === '/api/submit_query?q=B') return jsonResponse('task-B');
      if (url === '/api/fetch_query_result?task_id=task-A1') {
        return fetchMock.mock.calls.filter(([called]) => called === url).length === 1
          ? jsonResponse([], 202)
          : jsonResponse({}, 500);
      }
      const id = url.endsWith('task-A2') ? 'A' : 'B';
      return jsonResponse([{ reaction_id: id, proto: encodedReaction(id) }]);
    });
    vi.stubGlobal('fetch', fetchMock);
    const { result, rerender } = renderWithSharedClient('?q=A');
    await waitFor(() => expect(result.current.isError).toBe(true), { timeout: 3000 });
    rerender({ query: '?q=B' });
    await waitFor(() => expect(reactionIds(result.current.data)).toEqual(['B']));

    rerender({ query: '?q=A' });

    await waitFor(() => expect(reactionIds(result.current.data)).toEqual(['A']));
    expect(submitCalls(fetchMock)).toEqual([
      '/api/submit_query?q=A',
      '/api/submit_query?q=B',
      '/api/submit_query?q=A',
    ]);
  });

  describe('when a superseded search submit returns late', () => {
    it('polls its own task, not the current search task', async () => {
      const { fetchMock } = await raceSearches();
      expect(requestedUrls(fetchMock)).toContain(
        '/api/fetch_query_result?task_id=task-A',
      );
    });

    it('caches its own results for a return to that search', async () => {
      const { result, rerender } = await raceSearches();
      rerender({ query: '?q=A' });
      expect(reactionIds(result.current.data)).toEqual(['A']);
    });

    it('leaves the current search submitted once', async () => {
      const { fetchMock, result } = await raceSearches();
      await waitFor(() => expect(result.current.data?.status).toBe('success'), {
        timeout: 5000,
      });
      expect(reactionIds(result.current.data)).toEqual(['B']);
      expect(submitCalls(fetchMock)).toEqual([
        '/api/submit_query?q=A',
        '/api/submit_query?q=B',
      ]);
    });
  });
});
