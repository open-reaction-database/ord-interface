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
import reaction_pb from 'ord-schema';
import type { CompoundIdentifier } from 'ord-schema/proto/reaction_pb';
import { amountObj, amountStr } from '../../utils/amount';
import { enumName } from '../../utils/enum';
import { formattedTime } from '../../utils/outcomes';
import type { ReactionWorkupData } from '../../types/search';
import './WorkupsView.scss';

interface WorkupsViewProps {
  workup: ReactionWorkupData | undefined;
}

const getNameIdentifier = (identifiers: CompoundIdentifier.AsObject[]): string => {
  // Type 6 is NAME in CompoundIdentifierTypeMap; falls back to whatever
  // identifier is first if the compound wasn't given a name.
  const nameIdentifier = identifiers.find(identifier => identifier.type === 6);
  return nameIdentifier?.value ?? identifiers[0]?.value ?? '';
};

const WorkupsView: React.FC<WorkupsViewProps> = ({ workup }) => {
  if (!workup) return null;

  const workupType = enumName(reaction_pb.ReactionWorkup.ReactionWorkupType, workup.type) ?? '';
  const duration = formattedTime(workup.duration);
  const aliquotAmount = workup.amount ? amountStr(amountObj(workup.amount)) : '';

  return (
    <div className="workups-view">
      <div className="details">
        <div className="label">Type</div>
        <div className="value">{workupType}</div>

        {workup.details && (
          <>
            <div className="label">Details</div>
            <div className="value">{workup.details}</div>
          </>
        )}

        {duration && (
          <>
            <div className="label">Duration</div>
            <div className="value">{duration}</div>
          </>
        )}

        {aliquotAmount && (
          <>
            <div className="label">Aliquot amount</div>
            <div className="value">{aliquotAmount}</div>
          </>
        )}

        {workup.keepPhase && (
          <>
            <div className="label">Phase kept</div>
            <div className="value">{workup.keepPhase}</div>
          </>
        )}

        {workup.targetPh !== undefined && workup.targetPh !== 0 && (
          <>
            <div className="label">Target pH</div>
            <div className="value">{workup.targetPh}</div>
          </>
        )}

        {workup.isAutomated && (
          <>
            <div className="label">Automated</div>
            <div className="value">yes</div>
          </>
        )}
      </div>

      {workup.input && workup.input.componentsList.length > 0 && (
        <div className="inputs">
          <div className="title">Inputs</div>
          <div className="components">
            {workup.input.componentsList.map((component, idx) => (
              <React.Fragment key={idx}>
                <div className="identifier">{getNameIdentifier(component.identifiersList)}</div>
                <div className="amount">{amountStr(amountObj(component.amount))}</div>
              </React.Fragment>
            ))}
          </div>
        </div>
      )}
    </div>
  );
};

export default WorkupsView;
