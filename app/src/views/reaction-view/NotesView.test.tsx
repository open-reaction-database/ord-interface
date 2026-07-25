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

import { render } from '@testing-library/react';
import { describe, expect, it } from 'vitest';
import NotesView from './NotesView';
import type { ReactionNotesData } from '../../types/search';

const renderNotes = (notes: unknown): Record<string, string> => {
  const { container } = render(<NotesView notes={notes as ReactionNotesData} />);
  const labels = [...container.querySelectorAll('.label')].map(
    el => el.textContent ?? '',
  );
  const values = [...container.querySelectorAll('.value')].map(
    el => el.textContent ?? '',
  );
  return Object.fromEntries(labels.map((label, index) => [label, values[index] ?? '']));
};

describe('NotesView', () => {
  it('renders nothing without notes', () => {
    expect(renderNotes(undefined)).toEqual({});
  });

  // Field names come straight off the proto, so `isHeterogeneous` has to read
  // as "heterogeneous" and `procedureDetails` as "procedure details".
  it('turns camelCase field names into phrases', () => {
    expect(
      renderNotes({ isHeterogeneous: true, procedureDetails: 'stirred overnight' }),
    ).toEqual({
      heterogeneous: 'true',
      'procedure details': 'stirred overnight',
    });
  });

  // proto3 fills every unset scalar with false / "", which would otherwise
  // render a row per field in the message.
  it('skips fields left at their default', () => {
    expect(
      renderNotes({
        isHeterogeneous: false,
        formsPrecipitate: true,
        procedureDetails: '',
      }),
    ).toEqual({ 'forms precipitate': 'true' });
  });

  it('renders nothing when every field is at its default', () => {
    expect(renderNotes({ isHeterogeneous: false, procedureDetails: '' })).toEqual({});
  });
});
