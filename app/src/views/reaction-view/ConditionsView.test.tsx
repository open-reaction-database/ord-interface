/**
 * Copyright 2026 Open Reaction Database Project Authors
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

import { render } from '@testing-library/react';
import { describe, expect, it } from 'vitest';
import ConditionsView from './ConditionsView';
import type { ReactionConditionsData } from '../../types/search';

// Each pane is a grid of alternating .label/.value cells; read it back as a
// label-to-value map so the assertions read like the rendered table.
const fields = (container: HTMLElement): Record<string, string> => {
  const details = container.querySelector('.details');
  if (!details) return {};
  const labels = [...details.querySelectorAll('.label')].map(
    el => el.textContent ?? '',
  );
  const values = [...details.querySelectorAll('.value')].map(
    el => el.textContent ?? '',
  );
  return Object.fromEntries(labels.map((label, index) => [label, values[index] ?? '']));
};

const renderConditions = (conditions: unknown, display: string) =>
  fields(
    render(
      <ConditionsView
        conditions={conditions as ReactionConditionsData}
        display={display}
      />,
    ).container,
  );

describe('ConditionsView', () => {
  it('renders nothing without conditions', () => {
    const { container } = render(
      <ConditionsView
        conditions={undefined}
        display="temperature"
      />,
    );
    expect(container).toBeEmptyDOMElement();
  });

  it('renders nothing for an unknown pane', () => {
    const { container } = render(
      <ConditionsView
        conditions={{} as ReactionConditionsData}
        display="nonsense"
      />,
    );
    expect(container).toBeEmptyDOMElement();
  });

  describe('temperature', () => {
    it('shows the control type, details and setpoint', () => {
      expect(
        renderConditions(
          {
            temperature: {
              control: { type: 3, details: 'silicone oil' },
              setpoint: { value: 80, units: 1, precision: 0 },
            },
          },
          'temperature',
        ),
      ).toEqual({
        'Control Type': 'OIL_BATH',
        Details: 'silicone oil',
        Setpoint: '80 °C',
      });
    });

    it('reports a missing setpoint as "None"', () => {
      expect(renderConditions({ temperature: {} }, 'temperature')).toEqual({
        'Control Type': '',
        Setpoint: 'None',
      });
    });

    it('counts the recorded measurements', () => {
      expect(
        renderConditions({ temperature: { measurementsList: [{}, {}] } }, 'temperature')
          .Measurements,
      ).toBe('2 recorded');
    });
  });

  describe('pressure', () => {
    it('shows the control type, setpoint and atmosphere', () => {
      expect(
        renderConditions(
          {
            pressure: {
              control: { type: 5, details: 'balloon' },
              setpoint: { value: 3, units: 1, precision: 0 },
              atmosphere: { type: 6, details: 'from a balloon' },
            },
          },
          'pressure',
        ),
      ).toEqual({
        'Control Type': 'PRESSURIZED',
        Details: 'balloon',
        Setpoint: '3 bar',
        Atmosphere: 'HYDROGEN, from a balloon',
      });
    });

    it('counts the recorded measurements', () => {
      expect(
        renderConditions({ pressure: { measurementsList: [{}] } }, 'pressure')
          .Measurements,
      ).toBe('1 recorded');
    });
  });

  describe('stirring', () => {
    it('shows the type, rate and RPM', () => {
      expect(
        renderConditions(
          {
            stirring: {
              type: 3,
              details: 'egg-shaped bar',
              rate: { type: 1, rpm: 500 },
            },
          },
          'stirring',
        ),
      ).toEqual({
        Type: 'STIR_BAR',
        Details: 'egg-shaped bar',
        Rate: 'HIGH',
        RPM: '500',
      });
    });

    // An unset rate reads as UNSPECIFIED rather than as a blank cell.
    it('falls back to UNSPECIFIED without a rate', () => {
      expect(renderConditions({ stirring: { type: 3 } }, 'stirring').Rate).toBe(
        'UNSPECIFIED',
      );
    });

    // proto3 defaults rpm to 0, which is not a real stirring rate.
    it('hides a zero RPM', () => {
      expect(
        renderConditions(
          { stirring: { type: 3, rate: { type: 1, rpm: 0 } } },
          'stirring',
        ),
      ).not.toHaveProperty('RPM');
    });
  });

  describe('illumination', () => {
    it('shows the type, wavelength, color and distance', () => {
      expect(
        renderConditions(
          {
            illumination: {
              type: 4,
              details: '450 nm',
              color: 'blue',
              peakWavelength: { value: 450, units: 2, precision: 0 },
              distanceToVessel: { value: 5, units: 1, precision: 0 },
            },
          },
          'illumination',
        ),
      ).toEqual({
        Type: 'LED: 450 nm',
        'Peak Wavelength': '450 millimeter',
        Color: 'blue',
        'Distance to Vessel': '5 centimeter',
      });
    });

    it('reports missing lengths as "None"', () => {
      expect(renderConditions({ illumination: { type: 3 } }, 'illumination')).toEqual({
        Type: 'DARK',
        'Peak Wavelength': 'None',
        'Distance to Vessel': 'None',
      });
    });
  });

  describe('electrochemistry', () => {
    it('shows the type and the electrode materials', () => {
      expect(
        renderConditions(
          {
            electrochemistry: {
              type: 2,
              details: '10 mA',
              anodeMaterial: 'graphite',
              cathodeMaterial: 'nickel',
            },
          },
          'electrochemistry',
        ),
      ).toEqual({
        Type: 'CONSTANT_CURRENT',
        Details: '10 mA',
        Anode: 'graphite',
        Cathode: 'nickel',
      });
    });

    it('omits the electrodes when they were not recorded', () => {
      expect(
        renderConditions({ electrochemistry: { type: 3 } }, 'electrochemistry'),
      ).toEqual({ Type: 'CONSTANT_VOLTAGE' });
    });
  });

  describe('flow', () => {
    it('shows the type and pump', () => {
      expect(
        renderConditions(
          { flow: { type: 2, details: 'PFA tubing', pumpType: 'syringe' } },
          'flow',
        ),
      ).toEqual({
        Type: 'PLUG_FLOW_REACTOR',
        Details: 'PFA tubing',
        'Pump Type': 'syringe',
      });
    });
  });

  describe('other', () => {
    it('shows the flags and details that were set', () => {
      expect(
        renderConditions(
          {
            reflux: true,
            ph: 7.4,
            conditionsAreDynamic: true,
            details: 'ramped over 2 h',
          },
          'other',
        ),
      ).toEqual({
        Reflux: 'yes',
        pH: '7.4',
        'Conditions are dynamic': 'yes',
        Details: 'ramped over 2 h',
      });
    });

    // proto3 defaults pH to 0, which is indistinguishable from a strongly
    // acidic reading, so it is treated as unset.
    it('hides an unset pH', () => {
      expect(renderConditions({ ph: 0, reflux: true }, 'other')).toEqual({
        Reflux: 'yes',
      });
    });

    it('shows nothing when no other conditions were recorded', () => {
      expect(renderConditions({}, 'other')).toEqual({});
    });
  });
});
