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

import { render, screen } from '@testing-library/react';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import App from './App';

// Every route is reachable from a plain URL, so the router is exercised by
// pushing the path onto jsdom's history before rendering.
const renderAt = (path: string) => {
  window.history.pushState({}, '', path);
  return render(<App />);
};

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

  it('renders no page content for an unknown route', () => {
    const { container } = renderAt('/nope');
    expect(screen.getByRole('link', { name: 'Browse' })).toBeInTheDocument();
    expect(container.querySelector('.home')).toBeNull();
  });
});
