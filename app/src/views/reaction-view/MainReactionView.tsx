/**
 * Copyright 2023 Open Reaction Database Project Authors
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

import React, { useState, useEffect, useMemo, useCallback } from 'react';
import { useParams } from 'wouter';
import { Button, Paper, Title } from '@mantine/core';
import { ord } from 'ord-schema-protobufjs';
import clsx from 'clsx';
import CompoundView from './CompoundView';
import SetupView from './SetupView';
import ConditionsView from './ConditionsView';
import NotesView from './NotesView';
import ObservationsView from './ObservationsView';
import WorkupsView from './WorkupsView';
import OutcomesView from './OutcomesView';
import ProvenanceView from './ProvenanceView';
import EventsView from './EventsView';
import FloatingModal from '../../components/FloatingModal';
import LoadingSpinner from '../../components/LoadingSpinner';
import { api } from '../../utils/api';
import { decodeReaction, sortedInputEntries } from '../../utils/proto';
import { enumName } from '../../utils/enum';
import { formattedTime } from '../../utils/outcomes';
import type { ReactionData } from '../../types/search';
import '../../styles/tabs.scss';
import classes from './MainReactionView.module.scss';

const SETUP_TABS = ['vessel', 'environment', 'automation'];
const CONDITION_TABS = [
  'temperature',
  'pressure',
  'stirring',
  'illumination',
  'electrochemistry',
  'flow',
  'other',
];

const MainReactionView: React.FC = () => {
  const { reactionId } = useParams<{ reactionId: string }>();
  const [reaction, setReaction] = useState<ReactionData | null>(null);
  const [reactionSummary, setReactionSummary] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);
  const [inputsIdx, setInputsIdx] = useState(0);
  const [setupTab, setSetupTab] = useState('vessel');
  const [conditionTab, setConditionTab] = useState('temperature');
  const [workupsTab, setWorkupsTab] = useState(0);
  const [outcomesTab, setOutcomesTab] = useState(0);
  const [showRawReaction, setShowRawReaction] = useState(false);
  const [navItems, setNavItems] = useState<string[]>([]);
  const [activeNav, setActiveNav] = useState<string>('summary');

  const getReactionData = useCallback(async (): Promise<ReactionData | null> => {
    if (!reactionId) return null;

    try {
      const response = await api.post<Array<{ proto: string }>>('/reactions', {
        reaction_ids: [reactionId],
      });
      if (!response.data?.[0]?.proto) return null;
      return decodeReaction(response.data[0].proto);
    } catch (error) {
      console.error('Error fetching reaction data:', error);
      return null;
    }
  }, [reactionId]);

  const getReactionSummary = useCallback(async (): Promise<string> => {
    if (!reactionId) return '';

    try {
      const response = await api.get<string>('/reaction_summary', {
        params: { reaction_id: reactionId, compact: false },
      });
      return response.data;
    } catch (error) {
      // Don't pipe a 4xx/5xx HTML error body into the summary panel via
      // dangerouslySetInnerHTML; return empty so the section stays blank.
      console.error('Error fetching reaction summary:', error);
      return '';
    }
  }, [reactionId]);

  const inputEntries = useMemo(
    () => sortedInputEntries(reaction?.inputs),
    [reaction?.inputs],
  );

  const displayDetails = useMemo<Record<string, React.ReactNode>>(() => {
    const inputEntry = inputEntries[inputsIdx];
    if (!inputEntry) return {};
    const input: ord.IReactionInput = inputEntry[1];

    const raw = { ...(input as unknown as Record<string, unknown>) };
    const formatted: Record<string, React.ReactNode> = {};

    for (const [key, value] of Object.entries(raw)) {
      if (value == null || Array.isArray(value)) continue;
      formatted[key] = value as React.ReactNode;
    }

    if (input.additionDevice?.type != null) {
      formatted.additionDevice =
        enumName(
          ord.ReactionInput.AdditionDevice.AdditionDeviceType,
          input.additionDevice.type,
        )?.toLowerCase() ?? '';
    }

    if (input.additionSpeed?.type != null) {
      formatted.additionSpeed =
        enumName(
          ord.ReactionInput.AdditionSpeed.AdditionSpeedType,
          input.additionSpeed.type,
        )?.toLowerCase() ?? '';
    }

    if (input.additionDuration) {
      formatted.additionDuration = formattedTime(input.additionDuration) ?? '';
    }

    return formatted;
  }, [inputEntries, inputsIdx]);

  const displayConditionsOther = useMemo(() => {
    const conditions = reaction?.conditions;
    if (!conditions) return undefined;
    const otherFields: Array<keyof typeof conditions> = [
      'reflux',
      'ph',
      'conditionsAreDynamic',
      'details',
    ];
    return otherFields.find(key => conditions[key]);
  }, [reaction?.conditions]);

  const events = useMemo<ord.IRecordEvent[]>(() => {
    const provenance = reaction?.provenance;
    if (!provenance?.recordCreated) return [];

    const created: ord.IRecordEvent = {
      ...provenance.recordCreated,
      details: '(record created)',
    };
    const modified = provenance.recordModified ?? [];
    const all = [created, ...modified];

    all.sort((a, b) => {
      const aMs = a.time?.value ? new Date(a.time.value).getTime() : 0;
      const bMs = b.time?.value ? new Date(b.time.value).getTime() : 0;
      return aMs - bMs;
    });

    return all;
  }, [reaction?.provenance]);

  const getReactionType = (id: number | null | undefined): string =>
    enumName(ord.ReactionIdentifier.ReactionIdentifierType, id) ?? '';

  const getWorkupLabel = (type: number | null | undefined): string => {
    const name = enumName(ord.ReactionWorkup.ReactionWorkupType, type);
    return name ? name.toLowerCase().replace(/_/g, ' ') : '';
  };

  const scrollTo = (id: string) => {
    document.getElementById(id)?.scrollIntoView({ behavior: 'smooth' });
  };

  useEffect(() => {
    const fetchData = async () => {
      if (!reactionId) return;

      const [reactionData, summaryData] = await Promise.all([
        getReactionData(),
        getReactionSummary(),
      ]);

      if (reactionData) {
        setReaction(reactionData);
        setReactionSummary(summaryData);

        const items = ['summary', 'identifiers', 'inputs'];

        if (reactionData.setup) items.push('setup');
        if (reactionData.conditions) items.push('conditions');
        if (reactionData.notes) items.push('notes');
        if (reactionData.observations?.length) items.push('observations');
        if (reactionData.workups?.length) items.push('workups');

        items.push('outcomes', 'provenance', 'full-record');
        setNavItems(items);
      }

      setLoading(false);
    };

    fetchData();
  }, [reactionId, getReactionData, getReactionSummary]);

  useEffect(() => {
    const handleScroll = () => {
      const scrollPosition = window.scrollY + 120;

      for (let i = navItems.length - 1; i >= 0; i--) {
        const element = document.getElementById(navItems[i]);
        if (element) {
          const elementTop = element.offsetTop;
          if (scrollPosition >= elementTop) {
            setActiveNav(navItems[i]);
            break;
          }
        }
      }
    };

    window.addEventListener('scroll', handleScroll);
    handleScroll();
    return () => window.removeEventListener('scroll', handleScroll);
  }, [navItems]);

  if (loading) {
    return <LoadingSpinner />;
  }

  if (!reaction) {
    return (
      <Paper
        radius="sm"
        p="xl"
        className={classes.emptyState}
      >
        Failed to load reaction {reactionId}.
      </Paper>
    );
  }

  return (
    <div className={classes.reactionView}>
      <aside className={classes.navHolder}>
        <Paper
          radius="sm"
          p="sm"
          className={classes.nav}
        >
          {navItems.map(item => (
            <Button
              key={item}
              variant={activeNav === item ? 'filled' : 'transparent'}
              color={activeNav === item ? 'primary' : 'gray'}
              justify="flex-start"
              fullWidth
              size="compact-md"
              className={classes.navButton}
              onClick={() => scrollTo(item)}
            >
              {item.replace(/-/g, ' ')}
            </Button>
          ))}
        </Paper>
      </aside>

      <div className={classes.content}>
        <Paper
          id="summary"
          radius="sm"
          p="lg"
          className={classes.section}
        >
          <Title order={2}>Summary</Title>
          {reactionSummary && (
            <div className={classes.summaryScroll}>
              <div dangerouslySetInnerHTML={{ __html: reactionSummary }} />
            </div>
          )}
        </Paper>

        {reaction.identifiers.length > 0 && (
          <Paper
            id="identifiers"
            radius="sm"
            p="lg"
            className={classes.section}
          >
            <Title order={2}>Identifiers</Title>
            <div className={classes.identifiers}>
              {reaction.identifiers.map((identifier, index) => (
                <React.Fragment key={index}>
                  <div>{getReactionType(identifier.type)}</div>
                  <div>{identifier.value}</div>
                  <div>{identifier.details}</div>
                </React.Fragment>
              ))}
            </div>
          </Paper>
        )}

        {inputEntries.length > 0 && (
          <Paper
            id="inputs"
            radius="sm"
            p="lg"
            className={classes.section}
          >
            <Title order={2}>Inputs</Title>
            <div className="tabs">
              {inputEntries.map(([key], idx) => (
                <div
                  key={idx}
                  className={clsx('tab', { selected: inputsIdx === idx })}
                  onClick={() => setInputsIdx(idx)}
                >
                  {key}
                </div>
              ))}
            </div>
            <div className={classes.subTitle}>Details</div>
            <div className={classes.detailsGrid}>
              {Object.keys(displayDetails).map(key => (
                <React.Fragment key={key}>
                  <div className={classes.label}>{key.replace(/(?=[A-Z])/g, ' ')}</div>
                  <div>{displayDetails[key]}</div>
                </React.Fragment>
              ))}
            </div>
            <div className={classes.subTitle}>Components</div>
            <div className={classes.compounds}>
              {inputEntries[inputsIdx]?.[1].components?.map((component, index) => (
                <CompoundView
                  key={index}
                  component={component}
                />
              ))}
            </div>
          </Paper>
        )}

        {reaction.setup && (
          <Paper
            id="setup"
            radius="sm"
            p="lg"
            className={classes.section}
          >
            <Title order={2}>Setup</Title>
            <div className="tabs">
              {SETUP_TABS.map(
                tab =>
                  (tab !== 'automation' || reaction.setup?.isAutomated) && (
                    <div
                      key={tab}
                      className={clsx('tab', 'capitalize', {
                        selected: setupTab === tab,
                      })}
                      onClick={() => setSetupTab(tab)}
                    >
                      {tab}
                    </div>
                  ),
              )}
            </div>
            <SetupView
              setup={reaction.setup}
              display={setupTab}
            />
          </Paper>
        )}

        {reaction.conditions && (
          <Paper
            id="conditions"
            radius="sm"
            p="lg"
            className={classes.section}
          >
            <Title order={2}>Conditions</Title>
            <div className="tabs">
              {CONDITION_TABS.map(
                tab =>
                  (reaction.conditions?.[tab as keyof typeof reaction.conditions] ||
                    (tab === 'other' && displayConditionsOther)) && (
                    <div
                      key={tab}
                      className={clsx('tab', 'capitalize', {
                        selected: conditionTab === tab,
                      })}
                      onClick={() => setConditionTab(tab)}
                    >
                      {tab}
                    </div>
                  ),
              )}
            </div>
            <ConditionsView
              conditions={reaction.conditions}
              display={conditionTab}
            />
          </Paper>
        )}

        {reaction.notes && (
          <Paper
            id="notes"
            radius="sm"
            p="lg"
            className={classes.section}
          >
            <Title order={2}>Notes</Title>
            <NotesView notes={reaction.notes} />
          </Paper>
        )}

        {reaction.observations.length > 0 && (
          <Paper
            id="observations"
            radius="sm"
            p="lg"
            className={classes.section}
          >
            <Title order={2}>Observations</Title>
            <ObservationsView observations={reaction.observations} />
          </Paper>
        )}

        {reaction.workups.length > 0 && (
          <Paper
            id="workups"
            radius="sm"
            p="lg"
            className={classes.section}
          >
            <Title order={2}>Workups</Title>
            <div className="tabs">
              {reaction.workups.map((workup, idx) => (
                <div
                  key={idx}
                  className={clsx('tab', 'capitalize', {
                    selected: workupsTab === idx,
                  })}
                  onClick={() => setWorkupsTab(idx)}
                >
                  {getWorkupLabel(workup.type)}
                </div>
              ))}
            </div>
            <WorkupsView workup={reaction.workups[workupsTab]} />
          </Paper>
        )}

        {reaction.outcomes.length > 0 && (
          <Paper
            id="outcomes"
            radius="sm"
            p="lg"
            className={classes.section}
          >
            <Title order={2}>Outcomes</Title>
            <div className="tabs">
              {reaction.outcomes.map((_outcome, idx) => (
                <div
                  key={idx}
                  className={clsx('tab', 'capitalize', {
                    selected: outcomesTab === idx,
                  })}
                  onClick={() => setOutcomesTab(idx)}
                >
                  Outcome {idx + 1}
                </div>
              ))}
            </div>
            {/* key=outcomesTab so the inner tab/modal state resets when the user switches outcomes. */}
            <OutcomesView
              key={outcomesTab}
              outcome={reaction.outcomes[outcomesTab]}
            />
          </Paper>
        )}

        {reaction.provenance && (
          <Paper
            id="provenance"
            radius="sm"
            p="lg"
            className={classes.section}
          >
            <Title order={2}>Provenance</Title>
            <ProvenanceView provenance={reaction.provenance} />
          </Paper>
        )}

        {events.length > 0 && (
          <Paper
            id="events"
            radius="sm"
            p="lg"
            className={classes.section}
          >
            <Title order={2}>Record events</Title>
            <EventsView events={events} />
          </Paper>
        )}

        <Paper
          id="full-record"
          radius="sm"
          p="lg"
          className={classes.section}
        >
          <Title order={2}>Full record</Title>
          <div>
            <Button
              variant="default"
              radius="sm"
              onClick={() => setShowRawReaction(true)}
            >
              View full record
            </Button>
          </div>
        </Paper>

        {showRawReaction && (
          <FloatingModal
            title="Raw Data"
            onCloseModal={() => setShowRawReaction(false)}
          >
            <pre className={classes.rawRecord}>
              {JSON.stringify(ord.Reaction.toObject(reaction), null, 2)}
            </pre>
          </FloatingModal>
        )}
      </div>
    </div>
  );
};

export default MainReactionView;
