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

import { useQuery } from '@tanstack/react-query';
import { api } from '../utils/api';
import { decodeReaction } from '../utils/proto';
import type {
  NLInterpretation,
  NLQueryResponse,
  ResolvedComponent,
  SearchResult,
} from '../types/search';

export interface NLQueryData {
  interpretation: NLInterpretation;
  resolvedComponents: ResolvedComponent[];
  queryComponents: string[];
  results: SearchResult[];
  dryRun: boolean;
}

/**
 * Runs a natural-language search against `/api/nl_query`.
 *
 * Unlike structured search, the backend translates and executes in a single
 * synchronous request, so this is a plain fetch rather than the submit/poll
 * protocol of {@link useSearchTask}. Result protos are deserialized here the
 * same way, and the model's interpretation is surfaced so the page can show the
 * user how their question was understood.
 */
export function useNLQuery(query: string | null, enabled: boolean, dryRun = false) {
  return useQuery<NLQueryData>({
    queryKey: ['nl-query', query, dryRun],
    enabled: enabled && query !== null && query.trim() !== '',
    retry: false,
    staleTime: Infinity,
    queryFn: async (): Promise<NLQueryData> => {
      const { data: raw } = await api.get<NLQueryResponse>('/nl_query', {
        params: { q: query, ...(dryRun ? { dry_run: 'true' } : {}) },
      });
      const results: SearchResult[] = raw.results.map(r => ({
        ...r,
        data: decodeReaction(r.proto),
      }));
      return {
        interpretation: raw.interpretation,
        resolvedComponents: raw.resolved_components,
        queryComponents: raw.query_components,
        results,
        dryRun: raw.dry_run,
      };
    },
  });
}
