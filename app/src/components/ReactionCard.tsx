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
import './ReactionCard.scss';

interface ReactionCardProps {
  reaction: any;
  isSelectable?: boolean;
  isSelected?: boolean;
}

const ReactionCard: React.FC<ReactionCardProps> = ({ 
  reaction, 
  isSelectable = true, 
  isSelected = false 
}) => {
  return (
    <div className={`reaction-card ${isSelected ? 'reaction-card--selected' : ''}`}>
      <div className="reaction-card__header">
        <h3 className="reaction-card__title">
          Reaction ID: {reaction.reaction_id || 'Unknown'}
        </h3>
        {isSelectable && (
          <input 
            type="checkbox" 
            checked={isSelected}
            className="reaction-card__checkbox"
            readOnly
          />
        )}
      </div>
      <div className="reaction-card__content">
        <p className="reaction-card__description">
          {reaction.data?.identifiersList?.[0]?.value || 'No description available'}
        </p>
      </div>
    </div>
  );
};

export default ReactionCard;