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

/**
 * Fetches JSON from the API, throwing on a non-2xx response.
 *
 * Centralizes the "check response.ok, otherwise throw with the HTTP status"
 * boilerplate shared by the data-fetching hooks and components. The thrown
 * Error message is `${label} failed (HTTP ${status})`, where `label` defaults
 * to the request URL.
 *
 * Note: callers that intentionally *swallow* fetch errors — e.g. to keep an
 * HTML error body out of `dangerouslySetInnerHTML` — keep their own inline
 * `response.ok` handling and do not use this helper.
 */
export async function fetchJson<T>(
  input: RequestInfo | URL,
  init?: RequestInit,
  label?: string,
): Promise<T> {
  const response = await fetch(input, init);
  if (!response.ok) {
    const name = label ?? (typeof input === 'string' ? input : 'request');
    throw new Error(`${name} failed (HTTP ${response.status})`);
  }
  return (await response.json()) as T;
}
