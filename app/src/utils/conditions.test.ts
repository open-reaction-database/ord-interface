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

import type {
  IlluminationConditions,
  PressureConditions,
  StirringConditions,
} from 'ord-schema/proto/reaction_pb';
import { describe, expect, it } from 'vitest';
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
} from './conditions';

// The generated AsObject types narrow enum fields to their known values and
// require every scalar, so the fixtures below are cast in rather than typed.
// Several cases deliberately use a unit outside the enum.
const measurement = <T>(value: number, units: number, precision = 0): T =>
  ({ value, units, precision }) as T;

const enumValue = <T>(value: number): T => value as T;

const atmosphere = (
  type: number,
  details = '',
): PressureConditions.Atmosphere.AsObject =>
  ({ type, details }) as unknown as PressureConditions.Atmosphere.AsObject;

const stirringRate = (type: number): StirringConditions.StirringRate.AsObject =>
  ({
    type,
    details: '',
    rpm: 0,
  }) as unknown as StirringConditions.StirringRate.AsObject;

const illumination = (type: number, details = ''): IlluminationConditions.AsObject =>
  ({ type, details, color: '' }) as unknown as IlluminationConditions.AsObject;

describe('tempType', () => {
  it('names the control type', () => {
    expect(tempType(2)).toBe('AMBIENT');
    expect(tempType(6)).toBe('ICE_BATH');
  });

  it('renders an empty string when there is nothing to name', () => {
    expect(tempType(undefined)).toBe('');
    expect(tempType(99)).toBe('');
  });
});

describe('tempSetPoint', () => {
  it('renders "None" without a setpoint', () => {
    expect(tempSetPoint(undefined)).toBe('None');
  });

  // Only the unit's first letter is shown, after the degree sign: "°C".
  it('abbreviates the unit', () => {
    expect(tempSetPoint(measurement(25, 1))).toBe('25 °C');
    expect(tempSetPoint(measurement(77, 2))).toBe('77 °F');
    expect(tempSetPoint(measurement(298, 3))).toBe('298 °K');
  });

  it('includes the precision when one was recorded', () => {
    expect(tempSetPoint(measurement(25, 1, 2))).toBe('25 (± 2) °C');
  });

  it('omits the unit letter when the unit is unrecognized', () => {
    expect(tempSetPoint(measurement(25, 99))).toBe('25 °');
  });
});

describe('pressureType', () => {
  it('names the control type', () => {
    expect(pressureType(4)).toBe('SEALED');
  });

  it('renders an empty string when there is nothing to name', () => {
    expect(pressureType(undefined)).toBe('');
    expect(pressureType(99)).toBe('');
  });
});

describe('pressureSetPoint', () => {
  it('renders "None" without a setpoint', () => {
    expect(pressureSetPoint(undefined)).toBe('None');
  });

  it('spells the unit out in lowercase', () => {
    expect(pressureSetPoint(measurement(1, 2))).toBe('1 atmosphere');
    expect(pressureSetPoint(measurement(760, 8))).toBe('760 mm_hg');
  });

  it('includes the precision when one was recorded', () => {
    expect(pressureSetPoint(measurement(1, 1, 0.1))).toBe('1 (± 0.1) bar');
  });

  it('drops the unit when it is unrecognized', () => {
    expect(pressureSetPoint(measurement(1, 99))).toBe('1 ');
  });
});

describe('pressureAtmo', () => {
  it('names the atmosphere', () => {
    expect(pressureAtmo(atmosphere(3))).toBe('NITROGEN');
  });

  it('appends the details for a custom atmosphere', () => {
    expect(pressureAtmo(atmosphere(1, '5% H2 in Ar'))).toBe('CUSTOM, 5% H2 in Ar');
  });

  it('renders an empty string without an atmosphere', () => {
    expect(pressureAtmo(undefined)).toBe('');
  });
});

describe('stirType', () => {
  it('names the stirring method', () => {
    expect(stirType(3)).toBe('STIR_BAR');
  });

  it('renders an empty string when there is nothing to name', () => {
    expect(stirType(undefined)).toBe('');
    expect(stirType(99)).toBe('');
  });
});

describe('stirRate', () => {
  // The rate object is passed whole; the helper reads .type off it rather than
  // comparing the object itself to the enum values.
  it('names the rate from the rate object', () => {
    expect(stirRate(stirringRate(1))).toBe('HIGH');
    expect(stirRate(stirringRate(3))).toBe('LOW');
  });

  it('renders an empty string without a rate', () => {
    expect(stirRate(undefined)).toBe('');
  });
});

describe('illumType', () => {
  it('names the illumination type', () => {
    expect(illumType(illumination(4))).toBe('LED');
  });

  it('appends the details after a colon', () => {
    expect(illumType(illumination(4, '450 nm'))).toBe('LED: 450 nm');
  });

  it('renders an empty string without illumination', () => {
    expect(illumType(undefined)).toBe('');
  });
});

describe('lengthStr', () => {
  it('renders the value with a lowercase unit', () => {
    expect(lengthStr(measurement(5, 2))).toBe('5 millimeter');
  });

  it('drops the unit when it is unrecognized', () => {
    expect(lengthStr(measurement(5, 99))).toBe('5');
  });

  it('returns undefined without a length', () => {
    expect(lengthStr(undefined)).toBeUndefined();
  });
});

describe('electrochemType', () => {
  it('names the electrochemistry type', () => {
    expect(electrochemType(enumValue(2))).toBe('CONSTANT_CURRENT');
  });

  it('renders an empty string when there is nothing to name', () => {
    expect(electrochemType(undefined)).toBe('');
    expect(electrochemType(enumValue(99))).toBe('');
  });
});

describe('flowType', () => {
  it('names the flow type', () => {
    expect(flowType(enumValue(2))).toBe('PLUG_FLOW_REACTOR');
  });

  it('renders an empty string when there is nothing to name', () => {
    expect(flowType(undefined)).toBe('');
    expect(flowType(enumValue(99))).toBe('');
  });
});
