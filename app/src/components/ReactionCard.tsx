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
import reaction_pb from 'ord-schema';
import type { CompoundIdentifier, ProductMeasurement } from 'ord-schema/proto/reaction_pb';
import LoadingSpinner from './LoadingSpinner';
import CopyButton from './CopyButton';
import { enumName } from '../utils/enum';
import { formatPercentage } from '../utils/outcomes';
import type { SearchResult, ReactionData } from '../types/search';
import './ReactionCard.scss';

interface ReactionCardProps {
  reaction: SearchResult & { dataset_id?: string };
  isSelectable?: boolean;
  isSelected?: boolean;
  onSelectionChange?: (reactionId: string, isSelected: boolean) => void;
}

const YIELD_MEASUREMENT_TYPE = 3; // ord-schema ProductMeasurementType.YIELD

const ReactionCard: React.FC<ReactionCardProps> = ({
  reaction,
  isSelectable = true,
  isSelected = false,
  onSelectionChange,
}) => {
  const navigate = useNavigate();
  const [reactionTable, setReactionTable] = useState<string | null>(null);

  const getReactionTable = useCallback(async () => {
    try {
      const response = await fetch(`/api/reaction_summary?reaction_id=${reaction.reaction_id}`);
      // Skip the 4xx/5xx body — it's an HTML error page that
      // dangerouslySetInnerHTML would render verbatim in every card.
      if (!response.ok) {
        console.error(`reaction_summary failed (HTTP ${response.status}) for ${reaction.reaction_id}`);
        return;
      }
      setReactionTable(await response.text());
    } catch (error) {
      console.error('Error fetching reaction table:', error);
    }
  }, [reaction.reaction_id]);

  const getYield = (measurements: ProductMeasurement.AsObject[] = []): string => {
    const yieldObj = measurements.find(m => m.type === YIELD_MEASUREMENT_TYPE);
    return yieldObj?.percentage ? formatPercentage(yieldObj.percentage) : 'Not listed';
  };

  const getConversion = (data: ReactionData | undefined): string => {
    const conversion = data?.outcomesList?.[0]?.conversion;
    if (!conversion) return 'Not listed';
    return formatPercentage(conversion);
  };

  const conditionsAndDuration = (data: ReactionData | undefined): string[] => {
    const details: string[] = [];
    if (!data) return details;

    const temp = data.conditions?.temperature?.setpoint;
    if (temp) {
      const units = enumName(reaction_pb.Temperature.TemperatureUnit, temp.units);
      details.push(`at ${temp.value}${units ? ` ${units.toLowerCase()}` : '°C'}`);
    }

    const pressure = data.conditions?.pressure?.setpoint;
    if (pressure) {
      const units = enumName(reaction_pb.Pressure.PressureUnit, pressure.units);
      details.push(`under ${pressure.value}${units ? ` ${units.toLowerCase()}` : ' atm'}`);
    }

    const reactionTime = data.outcomesList?.[0]?.reactionTime;
    if (reactionTime?.value) {
      const units = enumName(reaction_pb.Time.TimeUnit, reactionTime.units);
      details.push(`for ${reactionTime.value}${units ? ` ${units.toLowerCase()}` : 's'}`);
    }

    return details;
  };

  const productIdentifier = (identifier: CompoundIdentifier.AsObject): string => {
    const type = enumName(reaction_pb.CompoundIdentifier.CompoundIdentifierType, identifier.type);
    return `${type ?? ''}: ${identifier.value}`;
  };

  const handleCheckboxChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    if (onSelectionChange) {
      onSelectionChange(reaction.reaction_id, event.target.checked);
    }
  };

  const handleViewDetails = () => {
    navigate(`/id/${reaction.reaction_id}`);
  };

  useEffect(() => {
    getReactionTable();
  }, [getReactionTable]);

  const reactionData = reaction.data;
  const firstOutcome = reactionData?.outcomesList?.[0];
  const firstProduct = firstOutcome?.productsList?.[0];
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
          {reactionTable ? <div dangerouslySetInnerHTML={{ __html: reactionTable }} /> : <LoadingSpinner />}
        </div>

        <div className="info">
          <div className="col full">
            <button onClick={handleViewDetails}>View Full Details</button>
          </div>

          <div className="col">
            <div className="yield">Yield: {getYield(firstProduct?.measurementsList || [])}</div>
            <div className="conversion">Conversion: {getConversion(reactionData)}</div>
            <div className="conditions">
              Conditions: {conditionsAndDuration(reactionData).join('; ') || 'Not Listed'}
            </div>
            {firstProductIdentifier && (
              <div className="smile">
                <CopyButton textToCopy={firstProductIdentifier.value || ''} />
                <div className="value">Product {productIdentifier(firstProductIdentifier)}</div>
              </div>
            )}
          </div>

          <div className="col">
            <div className="creator">
              Uploaded by {provenance?.recordCreated?.person?.name || 'Unknown'},{' '}
              {provenance?.recordCreated?.person?.organization || 'Unknown'}
            </div>
            <div className="date">
              Uploaded on{' '}
              {provenance?.recordCreated?.time?.value
                ? new Date(provenance.recordCreated.time.value).toLocaleDateString()
                : 'Unknown'}
            </div>
            <div className="doi">DOI: {provenance?.doi || 'Not available'}</div>
            {provenance?.publicationUrl && (
              <div className="publication">
                <a
                  href={provenance.publicationUrl}
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  Publication URL
                </a>
              </div>
            )}
            <div className="dataset">
              Dataset:{' '}
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
