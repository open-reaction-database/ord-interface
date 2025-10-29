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

import React, { useState } from 'react';
import './EventsView.scss';

interface EventsViewProps {
  events: any[];
}

const EventsView: React.FC<EventsViewProps> = ({ events }) => {
  const [eventIdx, setEventIdx] = useState(0);

  if (!events?.length) return null;

  const currentEvent = events[eventIdx];

  return (
    <div className="events-view">
      <div className="tabs">
        {events.map((event, idx) => (
          <div
            key={idx}
            className={`tab ${eventIdx === idx ? 'selected' : ''}`}
            onClick={() => setEventIdx(idx)}
          >
            {new Date(event.time?.value).toLocaleString()}
          </div>
        ))}
      </div>
      
      <div className="details">
        {currentEvent?.details && (
          <>
            <div className="label">Details</div>
            <div className="value">{currentEvent.details}</div>
          </>
        )}
        
        {currentEvent?.person?.username && (
          <>
            <div className="label">Username</div>
            <div className="value">{currentEvent.person.username}</div>
          </>
        )}
        
        {currentEvent?.person?.name && (
          <>
            <div className="label">Name</div>
            <div className="value">{currentEvent.person.name}</div>
          </>
        )}
        
        {currentEvent?.person?.orcid && (
          <>
            <div className="label">ORCID</div>
            <div className="value">{currentEvent.person.orcid}</div>
          </>
        )}
        
        {currentEvent?.person?.organization && (
          <>
            <div className="label">Organization</div>
            <div className="value">{currentEvent.person.organization}</div>
          </>
        )}
        
        {currentEvent?.person?.email && (
          <>
            <div className="label">Email</div>
            <div className="value">{currentEvent.person.email}</div>
          </>
        )}
      </div>
    </div>
  );
};

export default EventsView;