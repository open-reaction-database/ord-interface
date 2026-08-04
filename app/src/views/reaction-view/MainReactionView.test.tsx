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
import reaction_pb from 'ord-schema';
import type { DateTime, Reaction, ReactionInput } from 'ord-schema/proto/reaction_pb';
import { MemoryRouter, Route, Routes } from 'react-router-dom';
import { afterEach, describe, expect, it, vi } from 'vitest';
import MainReactionView from './MainReactionView';

const encode = (reaction: Reaction): string =>
  btoa(String.fromCharCode(...reaction.serializeBinary()));

/** A reaction carrying only the sections each test needs. */
const buildReaction = (build: (reaction: Reaction) => void = () => {}): string => {
  const reaction = new reaction_pb.Reaction();
  reaction.setReactionId('ord-1');
  build(reaction);
  return encode(reaction);
};

const namedInput = (order: number, smiles?: string): ReactionInput => {
  const input = new reaction_pb.ReactionInput();
  input.setAdditionOrder(order);
  if (smiles) {
    const component = input.addComponents();
    const identifier = component.addIdentifiers();
    identifier.setType(reaction_pb.CompoundIdentifier.CompoundIdentifierType.SMILES);
    identifier.setValue(smiles);
  }
  return input;
};

const timestamp = (value: string): DateTime => {
  const time = new reaction_pb.DateTime();
  time.setValue(value);
  return time;
};

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
        proto: buildReaction(reaction => {
          reaction.getInputsMap().set('m1_m2', namedInput(1, 'CCO'));
          reaction.setNotes(new reaction_pb.ReactionNotes());
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
        proto: buildReaction(reaction => {
          reaction.setSetup(new reaction_pb.ReactionSetup());
          reaction.setConditions(new reaction_pb.ReactionConditions());
          reaction.addObservations().setComment('turned yellow');
          reaction
            .addWorkups()
            .setType(reaction_pb.ReactionWorkup.ReactionWorkupType.WASH);
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
        proto: buildReaction(reaction => {
          const identifier = reaction.addIdentifiers();
          identifier.setType(
            reaction_pb.ReactionIdentifier.ReactionIdentifierType.REACTION_SMILES,
          );
          identifier.setValue('CCO>>CC=O');
          identifier.setDetails('from the paper');
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
        proto: buildReaction(reaction => {
          reaction.getInputsMap().set('added_second', namedInput(2, 'CC=O'));
          reaction.getInputsMap().set('added_first', namedInput(1, 'CCO'));
        }),
      });
      const { container } = renderReaction();

      await screen.findByText('Inputs');
      expect(tabsIn(container, 'inputs')).toEqual(['added_first', 'added_second']);
    });

    it('shows the selected input and switches on click', async () => {
      const user = userEvent.setup();
      stubApi({
        proto: buildReaction(reaction => {
          reaction.getInputsMap().set('first', namedInput(1, 'CCO'));
          const second = namedInput(2, 'CC=O');
          second.setAdditionOrder(2);
          reaction.getInputsMap().set('second', second);
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
        proto: buildReaction(reaction => {
          const input = namedInput(1, 'CCO');
          const device = new reaction_pb.ReactionInput.AdditionDevice();
          device.setType(
            reaction_pb.ReactionInput.AdditionDevice.AdditionDeviceType.SYRINGE,
          );
          input.setAdditionDevice(device);
          const speed = new reaction_pb.ReactionInput.AdditionSpeed();
          speed.setType(
            reaction_pb.ReactionInput.AdditionSpeed.AdditionSpeedType.DROPWISE,
          );
          input.setAdditionSpeed(speed);
          const duration = new reaction_pb.Time();
          duration.setValue(30);
          duration.setUnits(reaction_pb.Time.TimeUnit.MINUTE);
          input.setAdditionDuration(duration);
          reaction.getInputsMap().set('m1', input);
        }),
      });
      renderReaction();

      await screen.findByText('Inputs');
      expect(screen.getByText('syringe')).toBeInTheDocument();
      expect(screen.getByText('dropwise')).toBeInTheDocument();
      expect(screen.getByText('30 minute(s)')).toBeInTheDocument();
    });
  });

  describe('setup', () => {
    it('offers the automation tab only for an automated setup', async () => {
      stubApi({
        proto: buildReaction(reaction => {
          const setup = new reaction_pb.ReactionSetup();
          setup.setIsAutomated(true);
          reaction.setSetup(setup);
        }),
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
        proto: buildReaction(reaction => {
          reaction.setSetup(new reaction_pb.ReactionSetup());
        }),
      });
      const { container } = renderReaction();

      await screen.findByText('Setup');
      expect(tabsIn(container, 'setup')).toEqual(['vessel', 'environment']);
    });
  });

  describe('conditions', () => {
    it('offers a tab per populated condition', async () => {
      stubApi({
        proto: buildReaction(reaction => {
          const conditions = new reaction_pb.ReactionConditions();
          conditions.setTemperature(new reaction_pb.TemperatureConditions());
          conditions.setStirring(new reaction_pb.StirringConditions());
          reaction.setConditions(conditions);
        }),
      });
      const { container } = renderReaction();

      await screen.findByText('Conditions');
      expect(tabsIn(container, 'conditions')).toEqual(['temperature', 'stirring']);
    });

    // "other" has no message of its own; it appears when any of the loose
    // condition fields is set.
    it('offers the other tab for a loose condition field', async () => {
      stubApi({
        proto: buildReaction(reaction => {
          const conditions = new reaction_pb.ReactionConditions();
          conditions.setReflux(true);
          reaction.setConditions(conditions);
        }),
      });
      const { container } = renderReaction();

      await screen.findByText('Conditions');
      expect(tabsIn(container, 'conditions')).toEqual(['other']);
    });

    it('switches condition panes on click', async () => {
      const user = userEvent.setup();
      stubApi({
        proto: buildReaction(reaction => {
          const conditions = new reaction_pb.ReactionConditions();
          conditions.setTemperature(new reaction_pb.TemperatureConditions());
          const stirring = new reaction_pb.StirringConditions();
          stirring.setType(reaction_pb.StirringConditions.StirringMethodType.STIR_BAR);
          conditions.setStirring(stirring);
          reaction.setConditions(conditions);
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
        proto: buildReaction(reaction => {
          reaction
            .addWorkups()
            .setType(reaction_pb.ReactionWorkup.ReactionWorkupType.DRY_IN_VACUUM);
          reaction
            .addWorkups()
            .setType(reaction_pb.ReactionWorkup.ReactionWorkupType.FILTRATION);
        }),
      });
      const { container } = renderReaction();

      await screen.findByText('Workups');
      expect(tabsIn(container, 'workups')).toEqual(['dry in vacuum', 'filtration']);
    });

    it('switches workups on click', async () => {
      const user = userEvent.setup();
      stubApi({
        proto: buildReaction(reaction => {
          reaction
            .addWorkups()
            .setType(reaction_pb.ReactionWorkup.ReactionWorkupType.EXTRACTION);
          const second = reaction.addWorkups();
          second.setType(reaction_pb.ReactionWorkup.ReactionWorkupType.FILTRATION);
          second.setDetails('through celite');
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
        proto: buildReaction(reaction => {
          reaction.addOutcomes();
          reaction.addOutcomes();
        }),
      });
      const { container } = renderReaction();

      await screen.findByText('Outcomes');
      expect(tabsIn(container, 'outcomes')).toEqual(['Outcome 1', 'Outcome 2']);
    });
  });

  describe('record events', () => {
    it('labels the creation event and orders events oldest first', async () => {
      stubApi({
        proto: buildReaction(reaction => {
          const provenance = new reaction_pb.ReactionProvenance();
          const createdEvent = new reaction_pb.RecordEvent();
          createdEvent.setTime(timestamp('2021-01-01T00:00:00Z'));
          provenance.setRecordCreated(createdEvent);
          const modified = provenance.addRecordModified();
          modified.setTime(timestamp('2020-01-01T00:00:00Z'));
          modified.setDetails('backdated edit');
          reaction.setProvenance(provenance);
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
        proto: buildReaction(reaction => {
          const provenance = new reaction_pb.ReactionProvenance();
          provenance.setCity('Cambridge');
          reaction.setProvenance(provenance);
        }),
      });
      renderReaction();

      await screen.findByText('Provenance');
      expect(screen.queryByText('Record Events')).not.toBeInTheDocument();
    });
  });

  describe('the full record', () => {
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
