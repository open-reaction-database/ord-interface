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
import './MainReactionView.scss';

// Inline utility function for now
const outcomesUtil = {
  formattedTime: (duration: any): string => {
    if (!duration?.value) return '';
    return `${duration.value} ${duration.units || 'seconds'}`;
  }
};

const MainReactionView: React.FC = () => {
  const { reactionId } = useParams<{ reactionId: string }>();
  const [reaction, setReaction] = useState<any>({});
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

  const getReactionData = useCallback(async (): Promise<any> => {
    if (!reactionId) return null;
    
    try {
      const response = await fetch('/api/reactions', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ "reaction_ids": [reactionId] })
      });
      
      if (response.ok) {
        const data = await response.json();
        if (data && data[0]?.proto) {
          const base64string = data[0].proto;
          const bytes = base64ToBytes(base64string);
          const reaction = reaction_pb.Reaction.deserializeBinary(new Uint8Array(bytes)).toObject();
          
          // sort inputs by addition order
          if (reaction.inputsMap) {
            reaction.inputsMap.sort((a: any, b: any) => a[1].additionOrder - b[1].additionOrder);
          }
          
          return reaction;
        }
      }
      return null;
    } catch (error) {
      console.error('Error fetching reaction data:', error);
      return null;
    }
  }, [reactionId]);

  const getReactionSummary = useCallback(async (): Promise<string> => {
    if (!reactionId) return '';
    
    try {
      const response = await fetch(`/api/reaction_summary?reaction_id=${reactionId}&compact=false`);
      return await response.text();
    } catch (error) {
      console.error('Error fetching reaction summary:', error);
      return '';
    }
  }, [reactionId]);

  const displayDetails = useMemo(() => {
    if (!reaction?.inputsMap?.[inputsIdx]) return {};
    
    let returnArr = { ...reaction.inputsMap[inputsIdx][1] };
    
    // filter out null/undefined values and arrays
    let formattedDetails = Object.fromEntries(
      Object.entries(returnArr).filter(([_, v]) => v != null && !Array.isArray(v))
    );
    
    // some details are objects that need to be broken down for display
    if (formattedDetails.additionDevice && typeof formattedDetails.additionDevice === 'object' && (formattedDetails.additionDevice as any).type) {
      const deviceTypes = reaction_pb.ReactionInput.AdditionDevice.AdditionDeviceType;
      const device = Object.keys(deviceTypes).find(key => (deviceTypes as any)[key] === (formattedDetails.additionDevice as any).type);
      formattedDetails.additionDevice = device?.toLowerCase() || '';
    }
    
    if (formattedDetails.additionSpeed && typeof formattedDetails.additionSpeed === 'object' && (formattedDetails.additionSpeed as any).type) {
      const speedTypes = reaction_pb.ReactionInput.AdditionSpeed.AdditionSpeedType;
      const speed = Object.keys(speedTypes).find(key => (speedTypes as any)[key] === (formattedDetails.additionSpeed as any).type);
      formattedDetails.additionSpeed = speed?.toLowerCase() || '';
    }
    
    if (formattedDetails.additionDuration) {
      formattedDetails.additionDuration = outcomesUtil.formattedTime(formattedDetails.additionDuration);
    }
    
    return formattedDetails;
  }, [reaction, inputsIdx]);

  const displayConditionsOther = useMemo(() => {
    const otherFields = ['reflux', 'ph', 'conditions_are_dynamic', 'details'];
    return otherFields.find(key => reaction?.conditions?.[key]);
  }, [reaction?.conditions]);

  const events = useMemo(() => {
    const eventArray: any[] = [];
    if (!reaction?.provenance?.recordCreated) return eventArray;
    
    // add events to array
    eventArray.push(reaction.provenance.recordCreated);
    eventArray[0].details = "(record created)";
    eventArray.push(...(reaction.provenance.recordModifiedList || []));
    
    // sort by date to be safe
    eventArray.sort((a, b) => {
      const dateA = new Date(a.time?.value);
      const dateB = new Date(b.time?.value);
      return dateA.getTime() - dateB.getTime();
    });
    
    return eventArray;
  }, [reaction?.provenance]);

  const getReactionType = (id: number): string => {
    const identifiers = reaction_pb.ReactionIdentifier.ReactionIdentifierType;
    return Object.keys(identifiers).find(key => (identifiers as any)[key] === id) || '';
  };

  const getWorkupLabel = (type: number): string => {
    const workupTypes = reaction_pb.ReactionWorkup.ReactionWorkupType;
    return Object.keys(workupTypes).find(key => (workupTypes as any)[key] === type)?.toLowerCase().replace(/_/g, ' ') || '';
  };

  const setNavItemsFunction = useCallback(() => {
    let items = ['summary', 'identifiers', 'inputs'];
    const optionals = ['setup', 'conditions', 'notes', 'observations', 'workups'];
    
    optionals.forEach(item => {
      if (reaction[item] || reaction[`${item}List`]?.length) {
        items.push(item);
      }
    });
    
    const lastItems = ['outcomes', 'provenance', 'full-record'];
    items.push(...lastItems);
    
    return items;
  }, [reaction]);

  const scrollTo = (id: string) => {
    document.getElementById(id)?.scrollIntoView({ behavior: 'smooth' });
  };

  useEffect(() => {
    const fetchData = async () => {
      if (!reactionId) return;
      
      const [reactionData, summaryData] = await Promise.all([
        getReactionData(),
        getReactionSummary()
      ]);
      
      if (reactionData) {
        setReaction(reactionData);
        setReactionSummary(summaryData);
        setNavItems(setNavItemsFunction());
      }
      
      setLoading(false);
    };
    
    fetchData();
  }, [reactionId, getReactionData, getReactionSummary, setNavItemsFunction]);

  // Scroll listener to highlight active nav item
  useEffect(() => {
    const handleScroll = () => {
      const scrollPosition = window.scrollY + 100; // Offset for better UX
      
      // Find which section is currently in view
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

    // Add scroll listener
    window.addEventListener('scroll', handleScroll);
    
    // Initial call to set active nav on load
    handleScroll();
    
    // Cleanup
    return () => {
      window.removeEventListener('scroll', handleScroll);
    };
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
            {navItems.map((item) => (
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
          <div id="summary" className="section">
            {reactionSummary && (
              <div className="summary">
                <div className="display" dangerouslySetInnerHTML={{ __html: reactionSummary }} />
              </div>
            )}
          </div>

          {reaction?.identifiersList?.length > 0 && (
            <>
              <div id="identifiers">
                <div className="title">Identifiers</div>
                <div className="section">
                  <div className="identifiers">
                    {reaction.identifiersList.map((identifier: any, index: number) => (
                      <React.Fragment key={index}>
                        <div className="value">{getReactionType(identifier.type)}</div>
                        <div className="value">{identifier.value}</div>
                        <div className="value">{identifier.details}</div>
                      </React.Fragment>
                    ))}
                  </div>
                </div>
              </div>
            </>
          )}

          {reaction?.inputsMap?.length > 0 && (
            <>
              <div id="inputs">
                <div className="title">Inputs</div>
                <div className="section">
                  <div className="tabs">
                    {reaction.inputsMap.map(([key]: [string], idx: number) => (
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
                      {Object.keys(displayDetails).map((key) => (
                        <React.Fragment key={key}>
                          <div className="label">{key.replace(/(?=[A-Z])/g, ' ')}</div>
                          <div className="value">{displayDetails[key] as React.ReactNode}</div>
                        </React.Fragment>
                      ))}
                    </div>
                    <div className="title">Components</div>
                    <div className="compound">
                      {reaction.inputsMap[inputsIdx][1].componentsList?.map((component: any, index: number) => (
                        <CompoundView key={index} component={component} />
                      ))}
                    </div>
                  </div>
                </div>
              </div>
            </>
          )}

          {reaction?.setup && (
            <>
              <div id="setup">
                <div className="title">Setup</div>
                <div className="section">
                  <div className="tabs">
                    {setupTabs.map((tab) => (
                      (tab !== 'automation' || reaction.setup.is_automated) && (
                        <div
                          key={tab}
                          className={`tab capitalize ${setupTab === tab ? 'selected' : ''}`}
                          onClick={() => setSetupTab(tab)}
                        >
                          {tab}
                        </div>
                      )
                    ))}
                  </div>
                  <div className="details">
                    <SetupView setup={reaction.setup} display={setupTab} />
                  </div>
                </div>
              </div>
            </>
          )}

          {reaction?.conditions && (
            <>
              <div id="conditions">
                <div className="title">Conditions</div>
                <div className="section">
                  <div className="tabs">
                    {conditionTabs.map((tab) => (
                      (reaction.conditions[tab] || (tab === 'other' && displayConditionsOther)) && (
                        <div
                          key={tab}
                          className={`tab capitalize ${conditionTab === tab ? 'selected' : ''}`}
                          onClick={() => setConditionTab(tab)}
                        >
                          {tab}
                        </div>
                      )
                    ))}
                  </div>
                  <div className="details">
                    <ConditionsView conditions={reaction.conditions} display={conditionTab} />
                  </div>
                </div>
              </div>
            </>
          )}

          {reaction?.notes && (
            <>
              <div id="notes">
                <div className="title">Notes</div>
                <div className="section">
                  <div className="details">
                    <NotesView notes={reaction.notes} />
                  </div>
                </div>
              </div>
            </>
          )}

          {reaction?.observationsList?.length > 0 && (
            <>
              <div id="observations">
                <div className="title">Observations</div>
                <div className="section">
                  <div className="details">
                    <ObservationsView observations={reaction.observationsList} />
                  </div>
                </div>
              </div>
            </>
          )}

          {reaction?.workupsList?.length > 0 && (
            <>
              <div id="workups">
                <div className="title">Workups</div>
                <div className="section">
                  <div className="tabs">
                    {reaction.workupsList.map((workup: any, idx: number) => (
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
            </>
          )}

          {reaction?.outcomesList?.length > 0 && (
            <>
              <div id="outcomes">
                <div className="title">Outcomes</div>
                <div className="section">
                  <div className="tabs">
                    {reaction.outcomesList.map((_: any, idx: number) => (
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
                    <OutcomesView outcome={reaction.outcomesList[outcomesTab]} />
                  </div>
                </div>
              </div>
            </>
          )}

          {reaction?.provenance && (
            <>
              <div id="provenance">
                <div className="title">Provenance</div>
                <div className="section">
                  <ProvenanceView provenance={reaction.provenance} />
                </div>
              </div>
            </>
          )}

          {events?.length > 0 && (
            <>
              <div id="events">
                <div className="title">Record Events</div>
                <div className="section">
                  <EventsView events={events} />
                </div>
              </div>
            </>
          )}

          {reaction && (
            <>
              <div id="full-record">
                <div className="title">Full Record</div>
                <div className="section">
                  <div className="full-record button" onClick={() => setShowRawReaction(true)}>
                    View Full Record
                  </div>
                </div>
              </div>
            </>
          )}

          {showRawReaction && (
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