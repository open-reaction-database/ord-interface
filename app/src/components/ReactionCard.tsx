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

import React, { useEffect, useState, useCallback } from 'react';
import { useNavigate } from 'react-router-dom';
import LoadingSpinner from './LoadingSpinner';
import CopyButton from './CopyButton';
import reaction_pb from 'ord-schema';
import './ReactionCard.scss';

interface ReactionCardProps {
  reaction: any;
  isSelectable?: boolean;
  isSelected?: boolean;
  onSelectionChange?: (reactionId: string, isSelected: boolean) => void;
}

const ReactionCard: React.FC<ReactionCardProps> = ({ 
  reaction, 
  isSelectable = true, 
  isSelected = false,
  onSelectionChange
}) => {
  const navigate = useNavigate();
  const [reactionTable, setReactionTable] = useState<string | null>(null);

  const getReactionTable = useCallback(async () => {
    try {
      const response = await fetch(`/api/reaction_summary?reaction_id=${reaction.reaction_id}`);
      const responseData = await response.text();
      setReactionTable(responseData);
    } catch (error) {
      console.error('Error fetching reaction table:', error);
    }
  }, [reaction.reaction_id]);

  const getYield = (measurements: any[] = []) => {
    const yieldObj = measurements.find(m => m.type === 3); // ord-schema type 3 == "YIELD"
    if (yieldObj?.percentage) {
      return `${yieldObj.percentage.value}%`;
    }
    return "Not listed";
  };

  const getConversion = (reactionData: any) => {
    if (!reactionData.outcomesList?.[0]?.conversion) return "Not listed";
    // TODO: decode conversion properly
    return "Not listed";
  };

  const conditionsAndDuration = (reactionData: any) => {
    const details: string[] = [];
    
    // get temp - simplified for now
    const temp = reactionData.conditions?.temperature?.setpoint;
    if (temp) {
      details.push(`at ${temp.value}${temp.units || '°C'}`);
    }

    // get Pressure - simplified for now
    const pressure = reactionData.conditions?.pressure?.setpoint;
    if (pressure) {
      details.push(`under ${pressure.value}${pressure.units || ' atm'}`);
    }

    // get duration - simplified for now
    const reactionTime = reactionData.outcomesList?.[0]?.reactionTime;
    if (reactionTime?.value) {
      details.push(`for ${reactionTime.value}${reactionTime.units || 's'}`);
    }

    return details;
  };

  const productIdentifier = (identifier: any) => {
    const identifierTypes = reaction_pb.CompoundIdentifier.CompoundIdentifierType;
    const identifierType = Object.keys(identifierTypes).find(
      key => (identifierTypes as any)[key] === identifier.type
    );
    return `${identifierType}: ${identifier.value}`;
  };

  const handleCheckboxChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    if (onSelectionChange) {
      onSelectionChange(reaction.reaction_id, event.target.checked);
    }
  };

  const handleViewDetails = () => {
    navigate(`/reaction/${reaction.reaction_id}`);
  };

  useEffect(() => {
    getReactionTable();
  }, [getReactionTable]);

  const reactionData = reaction.data;
  const firstProduct = reactionData?.outcomesList?.[0]?.productsList?.[0];
  const firstProductIdentifier = firstProduct?.identifiersList?.[0];
  const provenance = reactionData?.provenance;

  return (
    <div className="reaction-container">
      <div className={`row ${isSelected ? 'selected' : ''}`}>
        {isSelectable && (
          <div className="select">
            <input
              type="checkbox"
              id={`select_${reaction.reaction_id}`}
              value={reaction.reaction_id}
              checked={isSelected}
              onChange={handleCheckboxChange}
            />
            <label htmlFor={`select_${reaction.reaction_id}`}>Select reaction</label>
          </div>
        )}
        
        {provenance?.isMined && (
          <div className="is-mined">
            <div className="is-mined-badge">Mined</div>
          </div>
        )}

        <div className="reaction-table">
          {reactionTable ? (
            <div dangerouslySetInnerHTML={{ __html: reactionTable }} />
          ) : (
            <LoadingSpinner />
          )}
        </div>

        <div className="info">
          <div className="col full">
            <button onClick={handleViewDetails}>
              View Full Details
            </button>
          </div>
          
          <div className="col">
            <div className="yield">
              Yield: {getYield(firstProduct?.measurementsList)}
            </div>
            <div className="conversion">
              Conversion: {getConversion(reactionData)}
            </div>
            <div className="conditions">
              Conditions: {conditionsAndDuration(reactionData).join("; ") || "Not Listed"}
            </div>
            {firstProductIdentifier && (
              <div className="smile">
                <CopyButton textToCopy={firstProductIdentifier.value} />
                <div className="value">
                  Product {productIdentifier(firstProductIdentifier)}
                </div>
              </div>
            )}
          </div>

          <div className="col">
            <div className="creator">
              Uploaded by {provenance?.recordCreated?.person?.name}, {provenance?.recordCreated?.person?.organization}
            </div>
            <div className="date">
              Uploaded on {provenance?.recordCreated?.time?.value ? 
                new Date(provenance.recordCreated.time.value).toLocaleDateString() : 'Unknown'}
            </div>
            <div className="doi">
              DOI: {provenance?.doi || 'Not available'}
            </div>
            {provenance?.publicationUrl && (
              <div className="publication">
                <a href={provenance.publicationUrl} target="_blank" rel="noopener noreferrer">
                  Publication URL
                </a>
              </div>
            )}
            <div className="dataset">
              Dataset: {' '}
              <a 
                href={`/search?dataset_id=${reaction.dataset_id}&limit=100`}
                target="_blank" 
                rel="noopener noreferrer"
              >
                {reaction.dataset_id}
              </a>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default ReactionCard;