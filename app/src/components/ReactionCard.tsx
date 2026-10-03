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
import {
  CompoundIdentifier_CompoundIdentifierTypeSchema,
  Pressure_PressureUnitSchema,
  ProductMeasurement_ProductMeasurementType,
  Temperature_TemperatureUnitSchema,
  Time_TimeUnitSchema,
  type CompoundIdentifier,
  type ProductMeasurement,
  type Reaction,
} from '@buf/open-reaction-database_ord-schema.bufbuild_es/ord-schema/proto/reaction_pb';
import LoadingSpinner from './LoadingSpinner';
import CopyButton from './CopyButton';
import { enumName } from '../utils/enum';
import { formatPercentage } from '../utils/outcomes';
import type { SearchResult } from '../types/search';
import './ReactionCard.scss';

interface ReactionCardProps {
  reaction: SearchResult;
  isSelectable?: boolean;
  isSelected?: boolean;
  onSelectionChange?: (reactionId: string, isSelected: boolean) => void;
}

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
      const response = await fetch(
        `/api/reaction_summary?reaction_id=${reaction.reaction_id}`,
      );
      // Skip the 4xx/5xx body — it's an HTML error page that
      // dangerouslySetInnerHTML would render verbatim in every card.
      if (!response.ok) {
        console.error(
          `reaction_summary failed (HTTP ${response.status}) for ${reaction.reaction_id}`,
        );
        return;
      }
      setReactionTable(await response.text());
    } catch (error) {
      console.error('Error fetching reaction table:', error);
    }
  }, [reaction.reaction_id]);

  const getYield = (measurements: ProductMeasurement[] = []): string => {
    const yieldObj = measurements.find(
      m => m.type === ProductMeasurement_ProductMeasurementType.YIELD,
    );
    return yieldObj?.value.case === 'percentage'
      ? formatPercentage(yieldObj.value.value)
      : 'Not listed';
  };

  const getConversion = (data: Reaction | undefined): string => {
    const conversion = data?.outcomes[0]?.conversion;
    if (!conversion) return 'Not listed';
    return formatPercentage(conversion);
  };

  // An unset setpoint value reads as 0, the proto3 default.
  const conditionsAndDuration = (data: Reaction | undefined): string[] => {
    const details: string[] = [];
    if (!data) return details;

    const temp = data.conditions?.temperature?.setpoint;
    if (temp) {
      const units = enumName(Temperature_TemperatureUnitSchema, temp.units);
      details.push(`at ${temp.value ?? 0}${units ? ` ${units.toLowerCase()}` : '°C'}`);
    }

    const pressure = data.conditions?.pressure?.setpoint;
    if (pressure) {
      const units = enumName(Pressure_PressureUnitSchema, pressure.units);
      details.push(
        `under ${pressure.value ?? 0}${units ? ` ${units.toLowerCase()}` : ' atm'}`,
      );
    }

    const reactionTime = data.outcomes[0]?.reactionTime;
    if (reactionTime?.value) {
      const units = enumName(Time_TimeUnitSchema, reactionTime.units);
      details.push(
        `for ${reactionTime.value}${units ? ` ${units.toLowerCase()}` : 's'}`,
      );
    }

    return details;
  };

  const productIdentifier = (identifier: CompoundIdentifier): string => {
    const type = enumName(
      CompoundIdentifier_CompoundIdentifierTypeSchema,
      identifier.type,
    );
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
  const firstOutcome = reactionData?.outcomes[0];
  const firstProduct = firstOutcome?.products[0];
  const firstProductIdentifier = firstProduct?.identifiers[0];
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
            <button onClick={handleViewDetails}>View Full Details</button>
          </div>

          <div className="col">
            <div className="yield">
              Yield: {getYield(firstProduct?.measurements || [])}
            </div>
            <div className="conversion">Conversion: {getConversion(reactionData)}</div>
            <div className="conditions">
              Conditions:{' '}
              {conditionsAndDuration(reactionData).join('; ') || 'Not Listed'}
            </div>
            {firstProductIdentifier && (
              <div className="smile">
                <CopyButton textToCopy={firstProductIdentifier.value || ''} />
                <div className="value">
                  Product {productIdentifier(firstProductIdentifier)}
                </div>
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
