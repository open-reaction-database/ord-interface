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
import { MemoryRouter } from 'react-router-dom';
import { describe, expect, it } from 'vitest';
import HeaderNav from './HeaderNav';
import MainFooter from './MainFooter';

describe('HeaderNav', () => {
  const renderNav = () =>
    render(
      <MemoryRouter>
        <HeaderNav />
      </MemoryRouter>,
    );

  it('links to the in-app pages', () => {
    renderNav();
    expect(screen.getByRole('link', { name: 'Browse' })).toHaveAttribute(
      'href',
      '/browse',
    );
    expect(screen.getByRole('link', { name: 'Search' })).toHaveAttribute(
      'href',
      '/search',
    );
    expect(screen.getByRole('link', { name: 'About' })).toHaveAttribute(
      'href',
      '/about',
    );
  });

  it('links out to the editor and the docs', () => {
    renderNav();
    expect(screen.getByRole('link', { name: 'Contribute' })).toHaveAttribute(
      'href',
      'https://app.open-reaction-database.org',
    );
    expect(screen.getByRole('link', { name: 'Docs' })).toHaveAttribute(
      'href',
      'https://docs.open-reaction-database.org',
    );
  });

  // /ask exists but is deliberately unlisted while it is in development.
  it('does not advertise the natural-language search', () => {
    renderNav();
    expect(screen.queryByRole('link', { name: /ask/i })).not.toBeInTheDocument();
  });
});

describe('MainFooter', () => {
  it('links to GitHub and the docs without leaking the referrer', () => {
    render(<MainFooter />);

    const github = screen.getByRole('link', { name: 'GitHub' });
    expect(github).toHaveAttribute('href', 'https://github.com/Open-Reaction-Database');
    expect(github).toHaveAttribute('rel', 'noopener noreferrer');
    expect(github).toHaveAttribute('target', '_blank');
    expect(screen.getByRole('link', { name: 'Documentation' })).toHaveAttribute(
      'href',
      'https://docs.open-reaction-database.org',
    );
  });
});
