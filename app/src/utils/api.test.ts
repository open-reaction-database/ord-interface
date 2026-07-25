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

import { afterEach, describe, expect, it, vi } from 'vitest';
import { fetchJson } from './api';

const mockFetch = (response: Partial<Response>) => {
  const fetchMock = vi.fn().mockResolvedValue(response as Response);
  vi.stubGlobal('fetch', fetchMock);
  return fetchMock;
};

afterEach(() => {
  vi.unstubAllGlobals();
});

describe('fetchJson', () => {
  it('returns the parsed body on success', async () => {
    mockFetch({ ok: true, status: 200, json: async () => ({ hits: 3 }) });
    await expect(fetchJson<{ hits: number }>('/api/thing')).resolves.toEqual({
      hits: 3,
    });
  });

  it('passes the request options through to fetch', async () => {
    const fetchMock = mockFetch({ ok: true, status: 200, json: async () => [] });
    const init = { method: 'POST', body: 'payload' };
    await fetchJson('/api/thing', init);
    expect(fetchMock).toHaveBeenCalledWith('/api/thing', init);
  });

  it('throws with the status on a non-2xx response', async () => {
    mockFetch({ ok: false, status: 500, json: async () => ({}) });
    await expect(fetchJson('/api/thing')).rejects.toThrow(
      '/api/thing failed (HTTP 500)',
    );
  });

  it('names the request by its label when one is given', async () => {
    mockFetch({ ok: false, status: 404, json: async () => ({}) });
    await expect(fetchJson('/api/thing', undefined, 'submit_query')).rejects.toThrow(
      'submit_query failed (HTTP 404)',
    );
  });

  // A URL or Request object has no useful default name, so it falls back to a
  // generic one rather than stringifying the input.
  it('falls back to "request" for a non-string input', async () => {
    mockFetch({ ok: false, status: 503, json: async () => ({}) });
    await expect(fetchJson(new URL('https://example.com/api'))).rejects.toThrow(
      'request failed (HTTP 503)',
    );
  });

  // A failed response body is typically an HTML error page; parsing it as JSON
  // would throw something unhelpful instead of the status.
  it('does not read the body of a failed response', async () => {
    const json = vi.fn();
    mockFetch({ ok: false, status: 500, json });
    await expect(fetchJson('/api/thing')).rejects.toThrow();
    expect(json).not.toHaveBeenCalled();
  });
});
