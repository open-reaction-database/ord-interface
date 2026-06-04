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

import React from 'react';
import { formattedTime } from '../../utils/outcomes';
import type { ReactionObservationData } from '../../types/search';
import './ObservationsView.scss';

interface ObservationsViewProps {
  observations: ReactionObservationData[];
}

const ObservationsView: React.FC<ObservationsViewProps> = ({ observations }) => {
  if (!observations.length) return null;

  return (
    <div className="observations-view">
      <div className="details">
        <div className="label">Time</div>
        <div className="label">Comment</div>
        {observations.map((obs, idx) => (
          <React.Fragment key={idx}>
            <div className="value">{formattedTime(obs.time) ?? ''}</div>
            <div className="value">{obs.comment}</div>
          </React.Fragment>
        ))}
      </div>
    </div>
  );
};

export default ObservationsView;
