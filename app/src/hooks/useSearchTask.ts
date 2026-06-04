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
import reaction_pb from 'ord-schema';
import { base64ToBytes } from '../utils/base64';
import type { SearchResult } from '../types/search';

const POLL_INTERVAL_MS = 1000;
const POLL_TIMEOUT_MS = 120_000;

type TaskState = { status: 'success'; results: SearchResult[] } | { status: 'pending'; taskId: string };

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
  const taskRef = useRef<TaskRef>({ queryString: null, taskId: null, submitPromise: null, startTime: 0 });

  return useQuery<TaskState>({
    queryKey: ['search-task', queryString],
    enabled: enabled && queryString !== null,
    retry: false,
    staleTime: Infinity,
    refetchInterval: query => (query.state.data?.status === 'pending' ? POLL_INTERVAL_MS : false),
    refetchIntervalInBackground: false,
    queryFn: async (): Promise<TaskState> => {
      if (!queryString) return { status: 'success', results: [] };

      // queryString changed since the last call — start fresh.
      if (taskRef.current.queryString !== queryString) {
        taskRef.current = { queryString, taskId: null, submitPromise: null, startTime: 0 };
      }

      if (taskRef.current.taskId === null) {
        if (!taskRef.current.submitPromise) {
          taskRef.current.startTime = Date.now();
          taskRef.current.submitPromise = (async () => {
            const submitRes = await fetch(`/api/submit_query${queryString}`);
            if (!submitRes.ok) {
              throw new Error(`submit_query failed (HTTP ${submitRes.status})`);
            }
            return (await submitRes.json()) as string;
          })();
        }
        try {
          taskRef.current.taskId = await taskRef.current.submitPromise;
        } finally {
          taskRef.current.submitPromise = null;
        }
      }

      if (Date.now() - taskRef.current.startTime > POLL_TIMEOUT_MS) {
        const id = taskRef.current.taskId;
        taskRef.current.taskId = null;
        throw new Error(`Search task ${id} timed out after ${POLL_TIMEOUT_MS / 1000}s`);
      }

      const res = await fetch(`/api/fetch_query_result?task_id=${taskRef.current.taskId}`);

      if (res.status === 200) {
        const raw = (await res.json()) as Omit<SearchResult, 'data'>[];
        const results: SearchResult[] = raw.map(r => ({
          ...r,
          data: reaction_pb.Reaction.deserializeBinary(new Uint8Array(base64ToBytes(r.proto))).toObject(),
        }));
        taskRef.current.taskId = null;
        return { status: 'success', results };
      }

      if (res.status === 202) {
        return { status: 'pending', taskId: taskRef.current.taskId };
      }

      const id = taskRef.current.taskId;
      taskRef.current.taskId = null;
      throw new Error(`Search task ${id} failed (HTTP ${res.status})`);
    },
  });
}
