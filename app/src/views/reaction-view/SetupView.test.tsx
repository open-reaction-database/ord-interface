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
import SetupView from './SetupView';
import type { ReactionSetupData } from '../../types/search';

const fields = (container: HTMLElement): Record<string, string> => {
  const details = container.querySelector('.details');
  if (!details) return {};
  const labels = [...details.querySelectorAll('.label')].map(
    el => el.textContent ?? '',
  );
  const values = [...details.querySelectorAll('.value')].map(
    el => el.textContent ?? '',
  );
  return Object.fromEntries(labels.map((label, index) => [label, values[index] ?? '']));
};

const renderSetup = (setup: unknown, display: string) =>
  fields(
    render(
      <SetupView
        setup={setup as ReactionSetupData}
        display={display}
      />,
    ).container,
  );

describe('SetupView', () => {
  it('renders no pane for an unknown display', () => {
    expect(renderSetup({}, 'nonsense')).toEqual({});
  });

  describe('vessel', () => {
    it('names the type, material and volume', () => {
      expect(
        renderSetup(
          {
            vessel: {
              type: 2,
              details: 'borosilicate',
              material: { type: 2 },
              volume: { value: 25, units: 2, precision: 0 },
            },
          },
          'vessel',
        ),
      ).toEqual({
        Type: 'ROUND_BOTTOM_FLASK',
        Details: 'borosilicate',
        Material: 'GLASS',
        Volume: '25 milliliter',
      });
    });

    it('joins the attachments and preparations', () => {
      expect(
        renderSetup(
          {
            vessel: {
              type: 2,
              attachmentsList: [
                { type: 3, details: 'rubber' },
                { type: 4, details: '' },
              ],
              preparationsList: [{ type: 4, details: '' }],
            },
          },
          'vessel',
        ),
      ).toMatchObject({
        Attachments: 'SEPTUM: rubber, CAP',
        Preparations: 'FLAME_DRIED',
      });
    });

    it('falls back to UNSPECIFIED for a missing material', () => {
      expect(renderSetup({ vessel: { type: 2 } }, 'vessel')).toEqual({
        Type: 'ROUND_BOTTOM_FLASK',
        Material: 'UNSPECIFIED',
        Volume: '',
      });
    });

    it('renders an empty pane without a vessel', () => {
      expect(renderSetup({}, 'vessel')).toEqual({
        Type: '',
        Material: 'UNSPECIFIED',
        Volume: '',
      });
    });
  });

  describe('environment', () => {
    it('names the environment type', () => {
      expect(
        renderSetup(
          { environment: { type: 2, details: 'in a glovebox' } },
          'environment',
        ),
      ).toEqual({ Type: 'FUME_HOOD', Details: 'in a glovebox' });
    });

    it('falls back to UNSPECIFIED when none was recorded', () => {
      expect(renderSetup({}, 'environment')).toEqual({ Type: 'UNSPECIFIED' });
    });
  });

  describe('automation', () => {
    it('names the platform', () => {
      expect(renderSetup({ automationPlatform: 'Chemspeed' }, 'automation')).toEqual({
        Platform: 'Chemspeed',
      });
    });

    it('falls back to UNSPECIFIED when none was recorded', () => {
      expect(renderSetup({}, 'automation')).toEqual({ Platform: 'UNSPECIFIED' });
    });
  });
});
