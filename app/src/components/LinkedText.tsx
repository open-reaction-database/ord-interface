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

import React, { type ReactNode } from 'react';
import './LinkedText.scss';

// Web addresses and bare DOIs ("10.1021/jacs.6c05959", with or without a "DOI:"
// prefix) as they appear in free text such as dataset descriptions. The URL
// branch comes first so that the DOI in https://doi.org/10.1021/... is consumed
// as part of the URL instead of matching on its own.
const LINK_PATTERN = /(?:https?:\/\/|www\.)[^\s<>"']+|\b10\.\d{4,9}\/[^\s<>"']+/;

const TRAILING_PUNCTUATION = /[.,;:!?'"]+$/;

const OPENING_BRACKETS: Record<string, string> = {
  ')': '(',
  ']': '[',
  '}': '{',
};

const countOccurrences = (text: string, character: string): number =>
  text.split(character).length - 1;

// Punctuation at the end of a match usually belongs to the surrounding prose
// ("...available at https://doi.org/10.1021/jacs.6c05959.") rather than to the
// link. A closing bracket is kept only when the link contains its opener, so
// that "(DOI: 10.15227/orgsyn.095.0080)" drops the paren but a Wikipedia-style
// ".../Foo_(bar)" keeps it.
const trimTrailingPunctuation = (link: string): string => {
  let trimmed = link.replace(TRAILING_PUNCTUATION, '');
  for (;;) {
    const closingBracket = trimmed.slice(-1);
    const openingBracket = OPENING_BRACKETS[closingBracket];
    if (
      !openingBracket ||
      countOccurrences(trimmed, openingBracket) >=
        countOccurrences(trimmed, closingBracket)
    ) {
      return trimmed;
    }
    trimmed = trimmed.slice(0, -1).replace(TRAILING_PUNCTUATION, '');
  }
};

const getHref = (link: string): string => {
  if (/^https?:\/\//i.test(link)) {
    return link;
  }
  if (/^www\./i.test(link)) {
    return `https://${link}`;
  }
  return `https://doi.org/${link}`;
};

interface LinkedTextProps {
  text?: string | null;
}

// Renders free text with any URLs and DOIs it contains as links.
const LinkedText: React.FC<LinkedTextProps> = ({ text }) => {
  if (!text) {
    return null;
  }
  const pattern = new RegExp(LINK_PATTERN, 'gi');
  const nodes: ReactNode[] = [];
  let cursor = 0;
  let match: RegExpExecArray | null;
  while ((match = pattern.exec(text)) !== null) {
    const link = trimTrailingPunctuation(match[0]);
    if (match.index > cursor) {
      nodes.push(text.slice(cursor, match.index));
    }
    nodes.push(
      <a
        key={match.index}
        className="linked-text"
        href={getHref(link)}
        target="_blank"
        rel="noopener noreferrer"
      >
        {link}
      </a>,
    );
    cursor = match.index + link.length;
    // Resume scanning at the punctuation left out of the link.
    pattern.lastIndex = cursor;
  }
  if (cursor < text.length) {
    nodes.push(text.slice(cursor));
  }
  return <>{nodes}</>;
};

export default LinkedText;
