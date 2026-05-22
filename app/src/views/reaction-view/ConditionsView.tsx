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
import {
  electrochemType,
  flowType,
  illumType,
  lengthStr,
  pressureAtmo,
  pressureSetPoint,
  pressureType,
  stirRate,
  stirType,
  tempSetPoint,
  tempType,
} from '../../utils/conditions';
import type { ReactionConditionsData } from '../../types/search';
import './ConditionsView.scss';

interface ConditionsViewProps {
  conditions: ReactionConditionsData | undefined;
  display: string;
}

const ConditionsView: React.FC<ConditionsViewProps> = ({ conditions, display }) => {
  if (!conditions) return null;

  if (display === 'temperature') {
    const t = conditions.temperature;
    return (
      <div className="conditions-view">
        <div className="temperature details">
          <div className="label">Control Type</div>
          <div className="value">{tempType(t?.control?.type)}</div>
          {t?.control?.details && (
            <>
              <div className="label">Details</div>
              <div className="value">{t.control.details}</div>
            </>
          )}
          <div className="label">Setpoint</div>
          <div className="value">{tempSetPoint(t?.setpoint)}</div>
          {t?.measurementsList && t.measurementsList.length > 0 && (
            <>
              <div className="label">Measurements</div>
              {/* TODO: render temperature measurements once the schema's shape is fully ported. */}
              <div className="value">{t.measurementsList.length} recorded</div>
            </>
          )}
        </div>
      </div>
    );
  }

  if (display === 'pressure') {
    const p = conditions.pressure;
    return (
      <div className="conditions-view">
        <div className="pressure details">
          <div className="label">Control Type</div>
          <div className="value">{pressureType(p?.control?.type)}</div>
          {p?.control?.details && (
            <>
              <div className="label">Details</div>
              <div className="value">{p.control.details}</div>
            </>
          )}
          <div className="label">Setpoint</div>
          <div className="value">{pressureSetPoint(p?.setpoint)}</div>
          <div className="label">Atmosphere</div>
          <div className="value">{pressureAtmo(p?.atmosphere)}</div>
          {p?.measurementsList && p.measurementsList.length > 0 && (
            <>
              <div className="label">Measurements</div>
              {/* TODO: render pressure measurements once the schema's shape is fully ported. */}
              <div className="value">{p.measurementsList.length} recorded</div>
            </>
          )}
        </div>
      </div>
    );
  }

  if (display === 'stirring') {
    const s = conditions.stirring;
    return (
      <div className="conditions-view">
        <div className="stirring details">
          <div className="label">Type</div>
          <div className="value">{stirType(s?.type)}</div>
          {s?.details && (
            <>
              <div className="label">Details</div>
              <div className="value">{s.details}</div>
            </>
          )}
          <div className="label">Rate</div>
          <div className="value">{stirRate(s?.rate) || 'UNSPECIFIED'}</div>
          {s?.rate?.rpm !== undefined && s.rate.rpm !== 0 && (
            <>
              <div className="label">RPM</div>
              <div className="value">{s.rate.rpm}</div>
            </>
          )}
        </div>
      </div>
    );
  }

  if (display === 'illumination') {
    const i = conditions.illumination;
    const distance = lengthStr(i?.distanceToVessel);
    return (
      <div className="conditions-view">
        <div className="illumination details">
          <div className="label">Type</div>
          <div className="value">{illumType(i)}</div>
          <div className="label">Peak Wavelength</div>
          <div className="value">{lengthStr(i?.peakWavelength) ?? 'None'}</div>
          {i?.color && (
            <>
              <div className="label">Color</div>
              <div className="value">{i.color}</div>
            </>
          )}
          <div className="label">Distance to Vessel</div>
          <div className="value">{distance ?? 'None'}</div>
        </div>
      </div>
    );
  }

  if (display === 'electrochemistry') {
    const e = conditions.electrochemistry;
    return (
      <div className="conditions-view">
        <div className="electro details">
          <div className="label">Type</div>
          <div className="value">{electrochemType(e?.type)}</div>
          {e?.details && (
            <>
              <div className="label">Details</div>
              <div className="value">{e.details}</div>
            </>
          )}
          {e?.anodeMaterial && (
            <>
              <div className="label">Anode</div>
              <div className="value">{e.anodeMaterial}</div>
            </>
          )}
          {e?.cathodeMaterial && (
            <>
              <div className="label">Cathode</div>
              <div className="value">{e.cathodeMaterial}</div>
            </>
          )}
          {/* TODO: render current, voltage, electrodeSeparation, and measurementsList. */}
        </div>
      </div>
    );
  }

  if (display === 'flow') {
    const f = conditions.flow;
    return (
      <div className="conditions-view">
        <div className="electro details">
          <div className="label">Type</div>
          <div className="value">{flowType(f?.type)}</div>
          {f?.details && (
            <>
              <div className="label">Details</div>
              <div className="value">{f.details}</div>
            </>
          )}
          {f?.pumpType && (
            <>
              <div className="label">Pump Type</div>
              <div className="value">{f.pumpType}</div>
            </>
          )}
          {/* TODO: render tubing dimensions if/when needed. */}
        </div>
      </div>
    );
  }

  if (display === 'other') {
    return (
      <div className="conditions-view">
        <div className="other details">
          {conditions.reflux && (
            <>
              <div className="label">Reflux</div>
              <div className="value">yes</div>
            </>
          )}
          {/* See WorkupsView for the proto3 zero-default tradeoff: 0 means
              "default unset" *or* a real, strongly-acidic pH. Match the Vue
              v-if='conditions.ph' behavior and hide on 0. */}
          {conditions.ph !== undefined && conditions.ph !== 0 && (
            <>
              <div className="label">pH</div>
              <div className="value">{conditions.ph}</div>
            </>
          )}
          {conditions.conditionsAreDynamic && (
            <>
              <div className="label">Conditions are dynamic</div>
              <div className="value">yes</div>
            </>
          )}
          {conditions.details && (
            <>
              <div className="label">Details</div>
              <div className="value">{conditions.details}</div>
            </>
          )}
        </div>
      </div>
    );
  }

  return null;
};

export default ConditionsView;
