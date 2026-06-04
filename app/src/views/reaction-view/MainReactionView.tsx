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
import { useParams } from 'react-router-dom';
import reaction_pb from 'ord-schema';
import type { ReactionInput, RecordEvent } from 'ord-schema/proto/reaction_pb';
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
import { base64ToBytes } from '../../utils/base64';
import { enumName } from '../../utils/enum';
import { formattedTime } from '../../utils/outcomes';
import type { ReactionData } from '../../types/search';
import './MainReactionView.scss';

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

  const setupTabs = ['vessel', 'environment', 'automation'];
  const conditionTabs = ['temperature', 'pressure', 'stirring', 'illumination', 'electrochemistry', 'flow', 'other'];

  const getReactionData = useCallback(async (): Promise<ReactionData | null> => {
    if (!reactionId) return null;

    try {
      const response = await fetch('/api/reactions', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ reaction_ids: [reactionId] }),
      });

      if (!response.ok) return null;
      const data = (await response.json()) as Array<{ proto: string }>;
      if (!data?.[0]?.proto) return null;

      const bytes = base64ToBytes(data[0].proto);
      const decoded = reaction_pb.Reaction.deserializeBinary(new Uint8Array(bytes)).toObject();

      // Sort inputs by their declared addition order.
      decoded.inputsMap?.sort((a, b) => (a[1].additionOrder ?? 0) - (b[1].additionOrder ?? 0));

      return decoded;
    } catch (error) {
      console.error('Error fetching reaction data:', error);
      return null;
    }
  }, [reactionId]);

  const getReactionSummary = useCallback(async (): Promise<string> => {
    if (!reactionId) return '';

    try {
      const response = await fetch(`/api/reaction_summary?reaction_id=${reactionId}&compact=false`);
      // Don't pipe a 4xx/5xx HTML error body into the summary panel via
      // dangerouslySetInnerHTML; return empty so the section stays blank.
      if (!response.ok) {
        console.error(`reaction_summary failed (HTTP ${response.status}) for ${reactionId}`);
        return '';
      }
      return await response.text();
    } catch (error) {
      console.error('Error fetching reaction summary:', error);
      return '';
    }
  }, [reactionId]);

  const displayDetails = useMemo<Record<string, React.ReactNode>>(() => {
    const inputEntry = reaction?.inputsMap?.[inputsIdx];
    if (!inputEntry) return {};
    const input: ReactionInput.AsObject = inputEntry[1];

    const raw = { ...(input as unknown as Record<string, unknown>) };
    const formatted: Record<string, React.ReactNode> = {};

    for (const [key, value] of Object.entries(raw)) {
      if (value == null || Array.isArray(value)) continue;
      formatted[key] = value as React.ReactNode;
    }

    if (input.additionDevice?.type !== undefined) {
      formatted.additionDevice =
        enumName(
          reaction_pb.ReactionInput.AdditionDevice.AdditionDeviceType,
          input.additionDevice.type,
        )?.toLowerCase() ?? '';
    }

    if (input.additionSpeed?.type !== undefined) {
      formatted.additionSpeed =
        enumName(reaction_pb.ReactionInput.AdditionSpeed.AdditionSpeedType, input.additionSpeed.type)?.toLowerCase() ??
        '';
    }

    if (input.additionDuration) {
      formatted.additionDuration = formattedTime(input.additionDuration) ?? '';
    }

    return formatted;
  }, [reaction, inputsIdx]);

  const displayConditionsOther = useMemo(() => {
    const conditions = reaction?.conditions;
    if (!conditions) return undefined;
    const otherFields: Array<keyof typeof conditions> = ['reflux', 'ph', 'conditionsAreDynamic', 'details'];
    return otherFields.find(key => conditions[key]);
  }, [reaction?.conditions]);

  const events = useMemo<RecordEvent.AsObject[]>(() => {
    const provenance = reaction?.provenance;
    if (!provenance?.recordCreated) return [];

    const created: RecordEvent.AsObject = { ...provenance.recordCreated, details: '(record created)' };
    const modified = provenance.recordModifiedList ?? [];
    const all = [created, ...modified];

    all.sort((a, b) => {
      const aMs = a.time?.value ? new Date(a.time.value).getTime() : 0;
      const bMs = b.time?.value ? new Date(b.time.value).getTime() : 0;
      return aMs - bMs;
    });

    return all;
  }, [reaction?.provenance]);

  const getReactionType = (id: number): string =>
    enumName(reaction_pb.ReactionIdentifier.ReactionIdentifierType, id) ?? '';

  const getWorkupLabel = (type: number): string => {
    const name = enumName(reaction_pb.ReactionWorkup.ReactionWorkupType, type);
    return name ? name.toLowerCase().replace(/_/g, ' ') : '';
  };

  const scrollTo = (id: string) => {
    document.getElementById(id)?.scrollIntoView({ behavior: 'smooth' });
  };

  useEffect(() => {
    const fetchData = async () => {
      if (!reactionId) return;

      const [reactionData, summaryData] = await Promise.all([getReactionData(), getReactionSummary()]);

      if (reactionData) {
        setReaction(reactionData);
        setReactionSummary(summaryData);

        const items = ['summary', 'identifiers', 'inputs'];
        const optionals: Array<keyof ReactionData> = [
          'setup',
          'conditions',
          'notes',
          'observationsList',
          'workupsList',
        ];

        optionals.forEach(key => {
          const value = reactionData[key];
          if (Array.isArray(value) ? value.length > 0 : Boolean(value)) {
            items.push(key.replace(/List$/, ''));
          }
        });

        items.push('outcomes', 'provenance', 'full-record');
        setNavItems(items);
      }

      setLoading(false);
    };

    fetchData();
  }, [reactionId, getReactionData, getReactionSummary]);

  useEffect(() => {
    const handleScroll = () => {
      const scrollPosition = window.scrollY + 100;

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
    return (
      <div className="main-reaction-view">
        <div className="loading">
          <LoadingSpinner />
        </div>
      </div>
    );
  }

  return (
    <div className="main-reaction-view">
      <div className="reaction-transition">
        <div className="nav-holder">
          <div className="nav">
            {navItems.map(item => (
              <div
                key={item}
                className={`nav-item ${activeNav === item ? 'active' : ''}`}
                onClick={() => scrollTo(item)}
              >
                {item.replace(/-/g, ' ')}
              </div>
            ))}
          </div>
        </div>

        <div className="content">
          <div className="title">Summary</div>
          <div
            id="summary"
            className="section"
          >
            {reactionSummary && (
              <div className="summary">
                <div
                  className="display"
                  dangerouslySetInnerHTML={{ __html: reactionSummary }}
                />
              </div>
            )}
          </div>

          {reaction?.identifiersList && reaction.identifiersList.length > 0 && (
            <div id="identifiers">
              <div className="title">Identifiers</div>
              <div className="section">
                <div className="identifiers">
                  {reaction.identifiersList.map((identifier, index) => (
                    <React.Fragment key={index}>
                      <div className="value">{getReactionType(identifier.type)}</div>
                      <div className="value">{identifier.value}</div>
                      <div className="value">{identifier.details}</div>
                    </React.Fragment>
                  ))}
                </div>
              </div>
            </div>
          )}

          {reaction?.inputsMap && reaction.inputsMap.length > 0 && (
            <div id="inputs">
              <div className="title">Inputs</div>
              <div className="section">
                <div className="tabs">
                  {reaction.inputsMap.map(([key], idx) => (
                    <div
                      key={idx}
                      className={`tab ${inputsIdx === idx ? 'selected' : ''}`}
                      onClick={() => setInputsIdx(idx)}
                    >
                      {key}
                    </div>
                  ))}
                </div>
                <div className="input">
                  <div className="title">Details</div>
                  <div className="details">
                    {Object.keys(displayDetails).map(key => (
                      <React.Fragment key={key}>
                        <div className="label">{key.replace(/(?=[A-Z])/g, ' ')}</div>
                        <div className="value">{displayDetails[key]}</div>
                      </React.Fragment>
                    ))}
                  </div>
                  <div className="title">Components</div>
                  <div className="compound">
                    {reaction.inputsMap[inputsIdx][1].componentsList?.map((component, index) => (
                      <CompoundView
                        key={index}
                        component={component}
                      />
                    ))}
                  </div>
                </div>
              </div>
            </div>
          )}

          {reaction?.setup && (
            <div id="setup">
              <div className="title">Setup</div>
              <div className="section">
                <div className="tabs">
                  {setupTabs.map(
                    tab =>
                      (tab !== 'automation' || reaction.setup?.isAutomated) && (
                        <div
                          key={tab}
                          className={`tab capitalize ${setupTab === tab ? 'selected' : ''}`}
                          onClick={() => setSetupTab(tab)}
                        >
                          {tab}
                        </div>
                      ),
                  )}
                </div>
                <div className="details">
                  <SetupView
                    setup={reaction.setup}
                    display={setupTab}
                  />
                </div>
              </div>
            </div>
          )}

          {reaction?.conditions && (
            <div id="conditions">
              <div className="title">Conditions</div>
              <div className="section">
                <div className="tabs">
                  {conditionTabs.map(
                    tab =>
                      (reaction.conditions?.[tab as keyof typeof reaction.conditions] ||
                        (tab === 'other' && displayConditionsOther)) && (
                        <div
                          key={tab}
                          className={`tab capitalize ${conditionTab === tab ? 'selected' : ''}`}
                          onClick={() => setConditionTab(tab)}
                        >
                          {tab}
                        </div>
                      ),
                  )}
                </div>
                <div className="details">
                  <ConditionsView
                    conditions={reaction.conditions}
                    display={conditionTab}
                  />
                </div>
              </div>
            </div>
          )}

          {reaction?.notes && (
            <div id="notes">
              <div className="title">Notes</div>
              <div className="section">
                <div className="details">
                  <NotesView notes={reaction.notes} />
                </div>
              </div>
            </div>
          )}

          {reaction?.observationsList && reaction.observationsList.length > 0 && (
            <div id="observations">
              <div className="title">Observations</div>
              <div className="section">
                <div className="details">
                  <ObservationsView observations={reaction.observationsList} />
                </div>
              </div>
            </div>
          )}

          {reaction?.workupsList && reaction.workupsList.length > 0 && (
            <div id="workups">
              <div className="title">Workups</div>
              <div className="section">
                <div className="tabs">
                  {reaction.workupsList.map((workup, idx) => (
                    <div
                      key={idx}
                      className={`tab capitalize ${workupsTab === idx ? 'selected' : ''}`}
                      onClick={() => setWorkupsTab(idx)}
                    >
                      {getWorkupLabel(workup.type)}
                    </div>
                  ))}
                </div>
                <div className="details">
                  <WorkupsView workup={reaction.workupsList[workupsTab]} />
                </div>
              </div>
            </div>
          )}

          {reaction?.outcomesList && reaction.outcomesList.length > 0 && (
            <div id="outcomes">
              <div className="title">Outcomes</div>
              <div className="section">
                <div className="tabs">
                  {reaction.outcomesList.map((_outcome, idx) => (
                    <div
                      key={idx}
                      className={`tab capitalize ${outcomesTab === idx ? 'selected' : ''}`}
                      onClick={() => setOutcomesTab(idx)}
                    >
                      Outcome {idx + 1}
                    </div>
                  ))}
                </div>
                <div className="details">
                  {/* key=outcomesTab so the inner tab/modal state resets when the user switches outcomes. */}
                  <OutcomesView
                    key={outcomesTab}
                    outcome={reaction.outcomesList[outcomesTab]}
                  />
                </div>
              </div>
            </div>
          )}

          {reaction?.provenance && (
            <div id="provenance">
              <div className="title">Provenance</div>
              <div className="section">
                <ProvenanceView provenance={reaction.provenance} />
              </div>
            </div>
          )}

          {events.length > 0 && (
            <div id="events">
              <div className="title">Record Events</div>
              <div className="section">
                <EventsView events={events} />
              </div>
            </div>
          )}

          {reaction && (
            <div id="full-record">
              <div className="title">Full Record</div>
              <div className="section">
                <div
                  className="full-record button"
                  onClick={() => setShowRawReaction(true)}
                >
                  View Full Record
                </div>
              </div>
            </div>
          )}

          {showRawReaction && reaction && (
            <FloatingModal
              title="Raw Data"
              onCloseModal={() => setShowRawReaction(false)}
            >
              <div className="data">
                <pre>{JSON.stringify(reaction, null, 2)}</pre>
              </div>
            </FloatingModal>
          )}
        </div>
      </div>
    </div>
  );
};

export default MainReactionView;
