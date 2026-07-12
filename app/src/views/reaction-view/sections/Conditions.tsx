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

// Ported from ord-app's ReactionView/Conditions, adapted to ord.I* messages
// (enum values resolve to their key strings, matching ord-app's converters).
import { Flex, Title } from '@mantine/core';
import { ord } from 'ord-schema-protobufjs';
import { RequiredOptionalFields } from '../../../components/display/RequiredOptionalFields';
import { enumName } from '../../../utils/enum';
import { reactionBoolean, vpu } from '../../../utils/formatting';
import { useReaction } from '../context';
import classes from './sections.module.scss';

export function Conditions() {
  const reaction = useReaction();
  const conditions: ord.IReactionConditions = reaction.conditions ?? {};

  return (
    <Flex direction="column">
      <Flex justify="space-between">
        <Flex
          align="center"
          gap="sm"
        >
          <Title order={2}>Conditions</Title>
        </Flex>
      </Flex>
      <Flex
        direction="column"
        gap="sm"
      >
        <RequiredOptionalFields
          entity={conditions}
          requiredFields={[
            { label: 'Reflux', render: c => reactionBoolean(c.reflux) },
            { label: 'pH', render: c => c.ph },
            {
              label: 'Dynamic conditions',
              render: c => reactionBoolean(c.conditionsAreDynamic),
            },
          ]}
        />
        <span className={classes.sectionLabel}>Temperature</span>
        <RequiredOptionalFields
          entity={conditions}
          requiredFields={[
            {
              label: 'Setpoint',
              render: ({ temperature }) =>
                temperature?.setpoint
                  ? vpu(temperature.setpoint, ord.Temperature.TemperatureUnit)
                  : 'UNSPECIFIED',
            },
          ]}
          optionalFields={[
            {
              label: 'Control',
              render: ({ temperature }) =>
                enumName(
                  ord.TemperatureConditions.TemperatureControl.TemperatureControlType,
                  temperature?.control?.type,
                ),
            },
            {
              label: 'Control Details',
              render: ({ temperature }) => temperature?.control?.details,
            },
          ]}
        />
        <span className={classes.sectionLabel}>Pressure</span>
        <RequiredOptionalFields
          entity={conditions}
          requiredFields={[
            {
              label: 'Setpoint',
              render: ({ pressure }) =>
                pressure?.setpoint
                  ? vpu(pressure.setpoint, ord.Pressure.PressureUnit)
                  : 'UNSPECIFIED',
            },
            {
              label: 'Control',
              render: ({ pressure }) =>
                enumName(
                  ord.PressureConditions.PressureControl.PressureControlType,
                  pressure?.control?.type,
                ) ?? 'UNSPECIFIED',
            },
            {
              label: 'Atmosphere',
              render: ({ pressure }) =>
                enumName(
                  ord.PressureConditions.Atmosphere.AtmosphereType,
                  pressure?.atmosphere?.type,
                ) ?? 'UNSPECIFIED',
            },
          ]}
        />
        <span className={classes.sectionLabel}>Stirring</span>
        <RequiredOptionalFields
          entity={conditions}
          requiredFields={[]}
          optionalFields={[
            {
              label: 'Method',
              render: ({ stirring }) =>
                enumName(ord.StirringConditions.StirringMethodType, stirring?.type),
            },
            {
              label: 'Details',
              render: ({ stirring }) => stirring?.details,
            },
            {
              label: 'Rate',
              render: ({ stirring }) =>
                enumName(
                  ord.StirringConditions.StirringRate.StirringRateType,
                  stirring?.rate?.type,
                ),
            },
            {
              label: 'Rate Details',
              render: ({ stirring }) => stirring?.rate?.details,
            },
            {
              label: 'RPM',
              render: ({ stirring }) => stirring?.rate?.rpm,
            },
          ]}
        />
        <span className={classes.sectionLabel}>Illumination</span>
        <RequiredOptionalFields
          entity={conditions}
          requiredFields={[
            {
              label: 'Type',
              render: ({ illumination }) =>
                enumName(
                  ord.IlluminationConditions.IlluminationType,
                  illumination?.type,
                ) ?? 'UNSPECIFIED',
            },
          ]}
          optionalFields={[]}
        />
        <span className={classes.sectionLabel}>Electrochemistry</span>
        <RequiredOptionalFields
          entity={conditions}
          requiredFields={[
            {
              label: 'Type',
              render: ({ electrochemistry }) =>
                enumName(
                  ord.ElectrochemistryConditions.ElectrochemistryType,
                  electrochemistry?.type,
                ) ?? 'UNSPECIFIED',
            },
          ]}
          optionalFields={[]}
        />
        <span className={classes.sectionLabel}>Flow</span>
        <RequiredOptionalFields
          entity={conditions}
          requiredFields={[
            {
              label: 'Type',
              render: ({ flow }) =>
                enumName(ord.FlowConditions.FlowType, flow?.type) ?? 'UNSPECIFIED',
            },
          ]}
          optionalFields={[]}
        />
      </Flex>
    </Flex>
  );
}
