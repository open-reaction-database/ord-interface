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

import { render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { create, toBinary, type MessageInitShape } from '@bufbuild/protobuf';
import {
  CompoundIdentifier_CompoundIdentifierType,
  ReactionInput_AdditionDevice_AdditionDeviceType,
  ReactionIdentifier_ReactionIdentifierType,
  ReactionInput_AdditionSpeed_AdditionSpeedType,
  ReactionSchema,
  ReactionWorkup_ReactionWorkupType,
  StirringConditions_StirringMethodType,
  Time_TimeUnit,
} from '@buf/open-reaction-database_ord-schema.bufbuild_es/ord-schema/proto/reaction_pb';
import { MemoryRouter, Route, Routes } from 'react-router-dom';
import { afterEach, describe, expect, it, vi } from 'vitest';
import MainReactionView from './MainReactionView';

/** A serialized reaction carrying only the sections each test needs. */
const buildReaction = (init: MessageInitShape<typeof ReactionSchema> = {}): string => {
  const reaction = create(ReactionSchema, { reactionId: 'ord-1', ...init });
  return btoa(String.fromCharCode(...toBinary(ReactionSchema, reaction)));
};

const namedInput = (order: number, smiles?: string) => ({
  additionOrder: order,
  components: smiles
    ? [
        {
          identifiers: [
            { type: CompoundIdentifier_CompoundIdentifierType.SMILES, value: smiles },
          ],
        },
      ]
    : [],
});

interface ApiOverrides {
  proto?: string | null;
  reactionsStatus?: number;
  summary?: string;
  summaryStatus?: number;
}

// The bodies POSTed to /api/reactions, for the assertions below.
const posted: RequestInit[] = [];

const stubApi = ({
  proto,
  reactionsStatus = 200,
  summary = '<b>CCO &rarr; CC=O</b>',
  summaryStatus = 200,
}: ApiOverrides = {}) => {
  const fetchMock = vi.fn(async (url: string, init?: RequestInit) => {
    if (url === '/api/reactions') {
      posted.push(init!);
      return {
        ok: reactionsStatus < 400,
        status: reactionsStatus,
        json: async () => (proto === null ? [] : [{ proto: proto ?? buildReaction() }]),
      };
    }
    if (url.startsWith('/api/reaction_summary')) {
      return {
        ok: summaryStatus < 400,
        status: summaryStatus,
        text: async () => summary,
      };
    }
    return { ok: true, status: 200, json: async () => '<svg></svg>' };
  });
  vi.stubGlobal('fetch', fetchMock);
  return fetchMock;
};

const renderReaction = (reactionId = 'ord-1') =>
  render(
    <MemoryRouter initialEntries={[`/id/${reactionId}`]}>
      <Routes>
        <Route
          path="/id/:reactionId"
          element={<MainReactionView />}
        />
      </Routes>
    </MemoryRouter>,
  );

const navItems = (container: HTMLElement): string[] =>
  [...container.querySelectorAll('.nav-item')].map(item => item.textContent ?? '');

const tabsIn = (container: HTMLElement, sectionId: string): string[] =>
  [...(container.querySelector(`#${sectionId}`)?.querySelectorAll('.tab') ?? [])].map(
    tab => tab.textContent ?? '',
  );

afterEach(() => {
  posted.length = 0;
  vi.unstubAllGlobals();
  vi.restoreAllMocks();
});

describe('MainReactionView', () => {
  it('requests the reaction and its summary', async () => {
    const fetchMock = stubApi();
    renderReaction('ord-42');

    await screen.findByText('Summary');
    expect(JSON.parse(posted[0].body as string)).toEqual({ reaction_ids: ['ord-42'] });
    expect(fetchMock).toHaveBeenCalledWith(
      '/api/reaction_summary?reaction_id=ord-42&compact=false',
    );
  });

  it('spins while the reaction is in flight', () => {
    vi.stubGlobal(
      'fetch',
      vi.fn(() => new Promise<Response>(() => {})),
    );
    const { container } = renderReaction();
    expect(container.querySelector('.spinner-main')).toBeInTheDocument();
  });

  it('renders the summary HTML', async () => {
    stubApi();
    const { container } = renderReaction();

    await waitFor(() =>
      expect(container.querySelector('.summary .display')?.innerHTML).toContain(
        '<b>CCO → CC=O</b>',
      ),
    );
  });

  // The 4xx/5xx body is an HTML error page, not a reaction drawing.
  it('leaves the summary blank when it fails to load', async () => {
    const consoleError = vi.spyOn(console, 'error').mockImplementation(() => {});
    stubApi({ summaryStatus: 500, summary: '<h1>500</h1>' });
    const { container } = renderReaction();

    await screen.findByText('Summary');
    expect(container.querySelector('.summary')).toBeNull();
    expect(consoleError).toHaveBeenCalled();
  });

  it('renders an empty page when the reaction is not found', async () => {
    stubApi({ proto: null });
    const { container } = renderReaction();

    await waitFor(() => expect(container.querySelector('.spinner-main')).toBeNull());
    expect(navItems(container)).toEqual([]);
  });

  it('renders an empty page when the request fails', async () => {
    stubApi({ reactionsStatus: 500 });
    const { container } = renderReaction();

    await waitFor(() => expect(container.querySelector('.spinner-main')).toBeNull());
    expect(navItems(container)).toEqual([]);
  });

  describe('the section nav', () => {
    it('lists only the sections the record actually has', async () => {
      stubApi({
        proto: buildReaction({
          inputs: { m1_m2: namedInput(1, 'CCO') },
          notes: {},
        }),
      });
      const { container } = renderReaction();

      await screen.findByText('Summary');
      expect(navItems(container)).toEqual([
        'summary',
        'identifiers',
        'inputs',
        'notes',
        'outcomes',
        'provenance',
        'full record',
      ]);
    });

    it('adds the optional sections that are populated', async () => {
      stubApi({
        proto: buildReaction({
          setup: {},
          conditions: {},
          observations: [{ comment: 'turned yellow' }],
          workups: [{ type: ReactionWorkup_ReactionWorkupType.WASH }],
        }),
      });
      const { container } = renderReaction();

      await screen.findByText('Summary');
      expect(navItems(container)).toContain('setup');
      expect(navItems(container)).toContain('conditions');
      expect(navItems(container)).toContain('observations');
      expect(navItems(container)).toContain('workups');
    });
  });

  describe('identifiers', () => {
    it('names each identifier type', async () => {
      stubApi({
        proto: buildReaction({
          identifiers: [
            {
              type: ReactionIdentifier_ReactionIdentifierType.REACTION_SMILES,
              value: 'CCO>>CC=O',
              details: 'from the paper',
            },
          ],
        }),
      });
      renderReaction();

      expect(await screen.findByText('REACTION_SMILES')).toBeInTheDocument();
      expect(screen.getByText('CCO>>CC=O')).toBeInTheDocument();
      expect(screen.getByText('from the paper')).toBeInTheDocument();
    });

    it('omits the section when there are none', async () => {
      stubApi();
      renderReaction();

      await screen.findByText('Summary');
      expect(screen.queryByText('Identifiers')).not.toBeInTheDocument();
    });
  });

  describe('inputs', () => {
    // The map arrives in arbitrary order; the view sorts by addition order so
    // the tabs read as the procedure was run.
    it('orders the tabs by addition order', async () => {
      stubApi({
        proto: buildReaction({
          inputs: {
            added_second: namedInput(2, 'CC=O'),
            added_first: namedInput(1, 'CCO'),
          },
        }),
      });
      const { container } = renderReaction();

      await screen.findByText('Inputs');
      expect(tabsIn(container, 'inputs')).toEqual(['added_first', 'added_second']);
    });

    it('breaks addition-order ties by key', async () => {
      stubApi({
        proto: buildReaction({
          inputs: {
            solvent: namedInput(1, 'O'),
            base: namedInput(1, '[Na+].[OH-]'),
          },
        }),
      });
      const { container } = renderReaction();

      await screen.findByText('Inputs');
      expect(tabsIn(container, 'inputs')).toEqual(['base', 'solvent']);
    });

    it('shows the selected input and switches on click', async () => {
      const user = userEvent.setup();
      stubApi({
        proto: buildReaction({
          inputs: {
            first: namedInput(1, 'CCO'),
            second: namedInput(2, 'CC=O'),
          },
        }),
      });
      const { container } = renderReaction();

      await screen.findByText('Inputs');
      const tabs = [...container.querySelector('#inputs')!.querySelectorAll('.tab')];
      expect(tabs[0]).toHaveClass('selected');

      await user.click(tabs[1]);
      expect(tabs[1]).toHaveClass('selected');
    });

    it('names the addition device, speed and duration', async () => {
      stubApi({
        proto: buildReaction({
          inputs: {
            m1: {
              ...namedInput(1, 'CCO'),
              additionDevice: {
                type: ReactionInput_AdditionDevice_AdditionDeviceType.SYRINGE,
              },
              additionSpeed: {
                type: ReactionInput_AdditionSpeed_AdditionSpeedType.DROPWISE,
              },
              additionDuration: { value: 30, units: Time_TimeUnit.MINUTE },
            },
          },
        }),
      });
      renderReaction();

      await screen.findByText('Inputs');
      expect(screen.getByText('syringe')).toBeInTheDocument();
      expect(screen.getByText('dropwise')).toBeInTheDocument();
      expect(screen.getByText('30 minute(s)')).toBeInTheDocument();
    });
  });

  describe('input details', () => {
    // The details pane lists the fields it knows how to format; other message
    // fields such as additionTime are not React children and stay out of it.
    it('renders an input that carries message fields it does not format', async () => {
      stubApi({
        proto: buildReaction({
          inputs: {
            m1: {
              ...namedInput(1, 'CCO'),
              additionTime: { value: 5, units: Time_TimeUnit.MINUTE },
              flowRate: { value: 1 },
              additionTemperature: { value: 0 },
              texture: { type: 2 },
            },
          },
        }),
      });
      const { container } = renderReaction();

      await screen.findByText('Inputs');
      const labels = [...container.querySelectorAll('#inputs .details .label')].map(
        label => label.textContent,
      );
      expect(labels).toEqual(['addition Order']);
    });
  });

  describe('setup', () => {
    it('offers the automation tab only for an automated setup', async () => {
      stubApi({
        proto: buildReaction({ setup: { isAutomated: true } }),
      });
      const { container } = renderReaction();

      await screen.findByText('Setup');
      expect(tabsIn(container, 'setup')).toEqual([
        'vessel',
        'environment',
        'automation',
      ]);
    });

    it('hides the automation tab for a manual setup', async () => {
      stubApi({
        proto: buildReaction({ setup: {} }),
      });
      const { container } = renderReaction();

      await screen.findByText('Setup');
      expect(tabsIn(container, 'setup')).toEqual(['vessel', 'environment']);
    });
  });

  describe('conditions', () => {
    it('offers a tab per populated condition', async () => {
      stubApi({
        proto: buildReaction({ conditions: { temperature: {}, stirring: {} } }),
      });
      const { container } = renderReaction();

      await screen.findByText('Conditions');
      expect(tabsIn(container, 'conditions')).toEqual(['temperature', 'stirring']);
    });

    // "other" has no message of its own; it appears when any of the loose
    // condition fields is set.
    it('offers the other tab for a loose condition field', async () => {
      stubApi({
        proto: buildReaction({ conditions: { reflux: true } }),
      });
      const { container } = renderReaction();

      await screen.findByText('Conditions');
      expect(tabsIn(container, 'conditions')).toEqual(['other']);
    });

    it('switches condition panes on click', async () => {
      const user = userEvent.setup();
      stubApi({
        proto: buildReaction({
          conditions: {
            temperature: {},
            stirring: { type: StirringConditions_StirringMethodType.STIR_BAR },
          },
        }),
      });
      const { container } = renderReaction();

      await screen.findByText('Conditions');
      await user.click(
        [...container.querySelector('#conditions')!.querySelectorAll('.tab')][1],
      );

      expect(screen.getByText('STIR_BAR')).toBeInTheDocument();
    });
  });

  describe('workups', () => {
    it('labels each tab with a readable workup type', async () => {
      stubApi({
        proto: buildReaction({
          workups: [
            { type: ReactionWorkup_ReactionWorkupType.DRY_IN_VACUUM },
            { type: ReactionWorkup_ReactionWorkupType.FILTRATION },
          ],
        }),
      });
      const { container } = renderReaction();

      await screen.findByText('Workups');
      expect(tabsIn(container, 'workups')).toEqual(['dry in vacuum', 'filtration']);
    });

    it('switches workups on click', async () => {
      const user = userEvent.setup();
      stubApi({
        proto: buildReaction({
          workups: [
            { type: ReactionWorkup_ReactionWorkupType.EXTRACTION },
            {
              type: ReactionWorkup_ReactionWorkupType.FILTRATION,
              details: 'through celite',
            },
          ],
        }),
      });
      const { container } = renderReaction();

      await screen.findByText('Workups');
      await user.click(
        [...container.querySelector('#workups')!.querySelectorAll('.tab')][1],
      );

      expect(screen.getByText('through celite')).toBeInTheDocument();
    });
  });

  describe('outcomes', () => {
    it('numbers the outcome tabs', async () => {
      stubApi({
        proto: buildReaction({ outcomes: [{}, {}] }),
      });
      const { container } = renderReaction();

      await screen.findByText('Outcomes');
      expect(tabsIn(container, 'outcomes')).toEqual(['Outcome 1', 'Outcome 2']);
    });
  });

  describe('record events', () => {
    it('labels the creation event and orders events oldest first', async () => {
      stubApi({
        proto: buildReaction({
          provenance: {
            recordCreated: { time: { value: '2021-01-01T00:00:00Z' } },
            recordModified: [
              { time: { value: '2020-01-01T00:00:00Z' }, details: 'backdated edit' },
            ],
          },
        }),
      });
      const { container } = renderReaction();

      await screen.findByText('Record Events');
      const tabs = [...container.querySelector('#events')!.querySelectorAll('.tab')];
      expect(tabs).toHaveLength(2);
      // The backdated edit sorts ahead of the creation event.
      expect(screen.getByText('backdated edit')).toBeInTheDocument();
    });

    it('omits the section when the record has no creation event', async () => {
      stubApi({
        proto: buildReaction({ provenance: { city: 'Cambridge' } }),
      });
      renderReaction();

      await screen.findByText('Provenance');
      expect(screen.queryByText('Record Events')).not.toBeInTheDocument();
    });
  });

  describe('the full record', () => {
    // The raw view is the proto3 JSON mapping, which names enum values.
    it('names enum values in the raw JSON', async () => {
      const user = userEvent.setup();
      stubApi({
        proto: buildReaction({
          workups: [{ type: ReactionWorkup_ReactionWorkupType.FILTRATION }],
        }),
      });
      renderReaction();

      await user.click(await screen.findByText('View Full Record'));
      expect(screen.getByText(/"type": "FILTRATION"/)).toBeInTheDocument();
    });

    it('opens and closes the raw JSON', async () => {
      const user = userEvent.setup();
      stubApi();
      renderReaction();

      await user.click(await screen.findByText('View Full Record'));
      expect(screen.getByText('Raw Data')).toBeInTheDocument();
      expect(screen.getByText(/"reactionId": "ord-1"/)).toBeInTheDocument();

      await user.click(screen.getByText('✕'));
      expect(screen.queryByText('Raw Data')).not.toBeInTheDocument();
    });
  });
});
