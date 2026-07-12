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

import axios from 'axios';

/**
 * Shared axios client for the ord-interface API (axios matches ord-app's HTTP
 * stack). In dev, Vite proxies `/api` to the FastAPI backend; in production,
 * nginx does the same.
 */
export const api = axios.create({ baseURL: '/api' });

/** Formats an unknown error as `${label} failed…` with HTTP status when available. */
export function apiErrorMessage(error: unknown, label: string): string {
  if (axios.isAxiosError(error) && error.response) {
    return `${label} failed (HTTP ${error.response.status})`;
  }
  return `${label} failed: ${error instanceof Error ? error.message : String(error)}`;
}
