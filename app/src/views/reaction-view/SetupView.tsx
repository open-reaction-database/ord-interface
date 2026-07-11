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

import React, { useMemo } from 'react';
import { ord } from 'ord-schema-protobufjs';
import { enumName } from '../../utils/enum';
import type { ReactionSetupData } from '../../types/search';
import './SetupView.scss';

interface SetupViewProps {
  setup: ReactionSetupData | undefined;
  display: string;
}

const SetupView: React.FC<SetupViewProps> = ({ setup, display }) => {
  const vessel = useMemo(() => setup?.vessel, [setup?.vessel]);

  const vesselType = useMemo(() => {
    if (!vessel?.type) return '';
    return enumName(ord.Vessel.VesselType, vessel.type) ?? '';
  }, [vessel?.type]);

  const vesselVolume = useMemo(() => {
    if (!vessel?.volume) return '';
    const label = enumName(ord.Volume.VolumeUnit, vessel.volume.units);
    return `${vessel.volume.value} ${label ? String(label).toLowerCase() : ''}`;
  }, [vessel?.volume]);

  const vesselAttachments = useMemo(() => {
    if (!vessel?.attachments?.length) return '';
    return vessel.attachments
      .map(attach => {
        const type = enumName(ord.VesselAttachment.VesselAttachmentType, attach.type);
        return `${String(type ?? '')}${attach.details ? `: ${attach.details}` : ''}`;
      })
      .join(', ');
  }, [vessel?.attachments]);

  const vesselPrep = useMemo(() => {
    if (!vessel?.preparations?.length) return '';
    return vessel.preparations
      .map(prep => {
        const type = enumName(ord.VesselPreparation.VesselPreparationType, prep.type);
        return `${String(type ?? '')}${prep.details ? `: ${prep.details}` : ''}`;
      })
      .join(', ');
  }, [vessel?.preparations]);

  const vesselMaterial = useMemo(() => {
    if (!vessel?.material?.type) return '';
    return String(
      enumName(ord.VesselMaterial.VesselMaterialType, vessel.material.type) ?? '',
    );
  }, [vessel?.material?.type]);

  const reactionEnv = useMemo(() => {
    const envVal = setup?.environment?.type;
    if (!envVal) return null;
    return (
      String(
        enumName(
          ord.ReactionSetup.ReactionEnvironment.ReactionEnvironmentType,
          envVal,
        ) ?? '',
      ) || null
    );
  }, [setup?.environment?.type]);

  return (
    <div className="setup-view">
      {display === 'vessel' && (
        <div className="vessel details">
          <div className="label">Type</div>
          <div className="value">{vesselType}</div>

          {vessel?.details && (
            <>
              <div className="label">Details</div>
              <div className="value">{vessel.details}</div>
            </>
          )}

          <div className="label">Material</div>
          <div className="value">{vesselMaterial || 'UNSPECIFIED'}</div>

          <div className="label">Volume</div>
          <div className="value">{vesselVolume}</div>

          {vessel?.attachments && vessel.attachments.length > 0 && (
            <>
              <div className="label">Attachments</div>
              <div className="value">{vesselAttachments}</div>
            </>
          )}

          {vessel?.preparations && vessel.preparations.length > 0 && (
            <>
              <div className="label">Preparations</div>
              <div className="value">{vesselPrep}</div>
            </>
          )}
        </div>
      )}

      {display === 'environment' && (
        <div className="environment details">
          <div className="label">Type</div>
          <div className="value">{reactionEnv || 'UNSPECIFIED'}</div>

          {setup?.environment?.details && (
            <>
              <div className="label">Details</div>
              <div className="value">{setup.environment.details}</div>
            </>
          )}
        </div>
      )}

      {display === 'automation' && (
        <div className="automated details">
          <div className="label">Platform</div>
          <div className="value">{setup?.automationPlatform || 'UNSPECIFIED'}</div>
        </div>
      )}
    </div>
  );
};

export default SetupView;
