/**
 * Copyright 2025 Open Reaction Database Project Authors
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

import { renderToStaticMarkup } from 'react-dom/server';
import { describe, expect, it } from 'vitest';
import LinkedText from './LinkedText';

const render = (text?: string | null): string =>
  renderToStaticMarkup(<LinkedText text={text} />);

const getLinks = (text: string): { text: string; href: string }[] =>
  [...render(text).matchAll(/<a\b[^>]*href="([^"]*)"[^>]*>([^<]*)<\/a>/g)].map(
    match => ({ href: match[1], text: match[2] }),
  );

// The text with its tags removed, for checking that nothing is lost or
// duplicated on the way through.
const getText = (text: string): string => render(text).replace(/<[^>]+>/g, '');

describe('LinkedText', () => {
  it('renders nothing without text', () => {
    expect(render(undefined)).toBe('');
    expect(render(null)).toBe('');
    expect(render('')).toBe('');
  });

  it('leaves text without links alone', () => {
    const text = 'Reactions from Figure 1 (Table 2), 2016,1, 658-666.';
    expect(render(text)).toBe(text);
  });

  it('links URLs', () => {
    expect(getLinks('Data from https://github.com/connorcoley/rexgen_direct')).toEqual([
      {
        text: 'https://github.com/connorcoley/rexgen_direct',
        href: 'https://github.com/connorcoley/rexgen_direct',
      },
    ]);
  });

  it('links bare hosts over https', () => {
    expect(getLinks('Performed by www.solvechemistry.com in London')).toEqual([
      { text: 'www.solvechemistry.com', href: 'https://www.solvechemistry.com' },
    ]);
  });

  it('links bare DOIs through doi.org, leaving any prefix as text', () => {
    expect(getLinks('Reactions from Figure 1 of DOI: 10.1021/jacs.8b01523')).toEqual([
      {
        text: '10.1021/jacs.8b01523',
        href: 'https://doi.org/10.1021/jacs.8b01523',
      },
    ]);
    expect(getLinks('Suzuki coupling from doi:10.1021/jacs.2c08592')).toEqual([
      {
        text: '10.1021/jacs.2c08592',
        href: 'https://doi.org/10.1021/jacs.2c08592',
      },
    ]);
  });

  it('treats a DOI inside a URL as part of that URL', () => {
    expect(getLinks('Available at https://doi.org/10.1021/jacs.6c05959')).toEqual([
      {
        text: 'https://doi.org/10.1021/jacs.6c05959',
        href: 'https://doi.org/10.1021/jacs.6c05959',
      },
    ]);
  });

  it('leaves trailing sentence punctuation out of the link', () => {
    expect(
      getLinks('The publication is available at https://doi.org/10.1021/jacs.6c05959.'),
    ).toEqual([
      {
        text: 'https://doi.org/10.1021/jacs.6c05959',
        href: 'https://doi.org/10.1021/jacs.6c05959',
      },
    ]);
    expect(getLinks('From 10.26434/chemrxiv-2024-22jrq, by the Doyle group')).toEqual([
      {
        text: '10.26434/chemrxiv-2024-22jrq',
        href: 'https://doi.org/10.26434/chemrxiv-2024-22jrq',
      },
    ]);
  });

  it('drops a closing bracket the link did not open', () => {
    expect(
      getLinks('Org. Synth. 2018, 95, 80-96 (DOI: 10.15227/orgsyn.095.0080)'),
    ).toEqual([
      {
        text: '10.15227/orgsyn.095.0080',
        href: 'https://doi.org/10.15227/orgsyn.095.0080',
      },
    ]);
  });

  it('keeps a closing bracket the link did open', () => {
    const url = 'https://en.wikipedia.org/wiki/Amine_(chemistry)';
    expect(getLinks(`See ${url}.`)).toEqual([{ text: url, href: url }]);
  });

  it('links every reference in the text', () => {
    const text =
      'Performed by SOLVE Chemistry (https://www.solvechemistry.com/) and published in https://doi.org/10.48550/arXiv.2506.07619.';
    expect(getLinks(text)).toEqual([
      {
        text: 'https://www.solvechemistry.com/',
        href: 'https://www.solvechemistry.com/',
      },
      {
        text: 'https://doi.org/10.48550/arXiv.2506.07619',
        href: 'https://doi.org/10.48550/arXiv.2506.07619',
      },
    ]);
    expect(getText(text)).toBe(text);
  });

  it('opens links in a new tab without leaking the referrer', () => {
    const html = render('See https://doi.org/10.1021/jacs.6c05959');
    expect(html).toContain('target="_blank"');
    expect(html).toContain('rel="noopener noreferrer"');
    expect(html).toContain('class="linked-text"');
  });
});
