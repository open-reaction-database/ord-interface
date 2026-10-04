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

import { useRef } from 'react';
import { useQuery } from '@tanstack/react-query';
import { fromBinary } from '@bufbuild/protobuf';
import { ReactionSchema } from '@buf/open-reaction-database_ord-schema.bufbuild_es/ord-schema/proto/reaction_pb';
import { base64ToBytes } from '../utils/base64';
import { fetchJson } from '../utils/api';
import type { SearchResult } from '../types/search';

const POLL_INTERVAL_MS = 1000;
const POLL_TIMEOUT_MS = 120_000;
// Bounds each request, so a stalled one fails the search instead of hanging it.
// Longer than the proxy's 60-second proxy_read_timeout (ord_interface/nginx.conf),
// so a slow response ends as the proxy's 504 rather than a client abort.
const REQUEST_TIMEOUT_MS = 90_000;

type TaskState =
  | { status: 'success'; results: SearchResult[] }
  | { status: 'pending'; taskId: string };

interface TaskRef {
  queryString: string | null;
  taskId: string | null;
  // In-flight submit_query promise. Two queryFn invocations can overlap
  // (StrictMode dev double-invoke, react-query refetch racing a polling
  // refetch, etc.); holding the in-flight promise on the ref lets the second
  // caller await the first one's result instead of firing a duplicate
  // submit_query that leaves an orphaned task on the backend.
  submitPromise: Promise<string> | null;
  startTime: number;
}

/** The React Query key under which `useSearchTask` caches a query string's search. */
export const searchTaskKey = (queryString: string | null) =>
  ['search-task', queryString] as const;

/**
 * Runs the API's submit-query / poll-result protocol against the given query
 * string, returning the materialized search results once the task completes.
 *
 * Submits to `/api/submit_query<queryString>` exactly once per `queryString`
 * change, then polls `/api/fetch_query_result?task_id=…` every second until the
 * server returns 200. Gives up after `POLL_TIMEOUT_MS` to bound user wait.
 */
export function useSearchTask(queryString: string | null, enabled: boolean) {
  // All polling state lives on this single ref, keyed by the queryString that
  // owns it, so a queryString change is the *only* reset signal. An earlier
  // version reset taskId / startTime from a separate useEffect; under
  // <StrictMode> dev the effect was re-running after queryFn had already set
  // startTime, leaving it at 0 — and "Date.now() - 0 > 120s" tripped the
  // timeout on the very first poll iteration.
  const taskRef = useRef<TaskRef>({
    queryString: null,
    taskId: null,
    submitPromise: null,
    startTime: 0,
  });

  return useQuery<TaskState>({
    queryKey: searchTaskKey(queryString),
    enabled: enabled && queryString !== null,
    retry: false,
    staleTime: Infinity,
    // A failed or timed-out poll keeps the last pending result as data, so the
    // error status is what stops the polling.
    refetchInterval: query =>
      query.state.status !== 'error' && query.state.data?.status === 'pending'
        ? POLL_INTERVAL_MS
        : false,
    refetchIntervalInBackground: false,
    queryFn: async (): Promise<TaskState> => {
      if (!queryString) return { status: 'success', results: [] };

      // queryString changed since the last call — start fresh.
      if (taskRef.current.queryString !== queryString) {
        taskRef.current = {
          queryString,
          taskId: null,
          submitPromise: null,
          startTime: 0,
        };
      }

      // This call's own state. After an await the ref may belong to a newer
      // queryString, and a late call must not poll or clear that query's task.
      const task = taskRef.current;

      if (task.taskId === null) {
        if (!task.submitPromise) {
          task.startTime = Date.now();
          task.submitPromise = fetchJson<string>(
            `/api/submit_query${queryString}`,
            { signal: AbortSignal.timeout(REQUEST_TIMEOUT_MS) },
            'submit_query',
          );
        }
        try {
          task.taskId = await task.submitPromise;
        } finally {
          task.submitPromise = null;
        }
      }

      try {
        const res = await fetch(`/api/fetch_query_result?task_id=${task.taskId}`, {
          signal: AbortSignal.timeout(REQUEST_TIMEOUT_MS),
        });

        if (res.status === 200) {
          const raw = (await res.json()) as Omit<SearchResult, 'data'>[];
          const results: SearchResult[] = raw.map(r => ({
            ...r,
            data: fromBinary(ReactionSchema, new Uint8Array(base64ToBytes(r.proto))),
          }));
          task.taskId = null;
          return { status: 'success', results };
        }

        // The deadline applies only to a task that is still running, so a result
        // that is ready by the first poll past it is returned, not discarded.
        if (res.status === 202) {
          if (Date.now() - task.startTime > POLL_TIMEOUT_MS) {
            throw new Error(
              `Search task ${task.taskId} timed out after ${POLL_TIMEOUT_MS / 1000}s`,
            );
          }
          return { status: 'pending', taskId: task.taskId };
        }

        throw new Error(`Search task ${task.taskId} failed (HTTP ${res.status})`);
      } catch (error) {
        // Running a failed search again submits it afresh. The task may have
        // expired, may never finish, or may take too long to read again.
        task.taskId = null;
        throw error;
      }
    },
  });
}
