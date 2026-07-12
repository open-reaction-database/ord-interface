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

// Ported from ord-app's ReactionView/Workups (read-only). ord-app gates fields
// by workup type via WorkupConstants; here optional fields simply render when
// present, which matches for real data.
import { Accordion, Flex, Title } from '@mantine/core';
import { ord } from 'ord-schema-protobufjs';
import { Counter } from '../../../components/display/Counter';
import { RequiredOptionalFields } from '../../../components/display/RequiredOptionalFields';
import { enumName } from '../../../utils/enum';
import { formatAmount, reactionBoolean, vpu } from '../../../utils/formatting';
import { useReaction } from '../context';
import { InputComponentsListItem, InputsHeaderRow } from './Inputs';

export function Workups() {
  const reaction = useReaction();
  const workups = reaction.workups ?? [];

  return (
    <Flex
      direction="column"
      gap="md"
    >
      <Flex justify="space-between">
        <Flex
          align="center"
          gap="sm"
        >
          <Title order={2}>Workups</Title>
          <Counter amount={workups.length} />
        </Flex>
      </Flex>
      <span>
        Workup steps refer to any additions, purifications, or other operations after
        the 'reaction' stage prior to analysis
      </span>
      <Flex
        direction="column"
        gap="sm"
      >
        {workups.map((workup, index) => (
          <Flex
            direction="column"
            key={index}
          >
            <Flex align="center">
              <Title order={3}>Workup {index + 1}</Title>
            </Flex>
            <RequiredOptionalFields
              entity={workup}
              requiredFields={[
                {
                  label: 'Type',
                  render: item =>
                    enumName(ord.ReactionWorkup.ReactionWorkupType, item.type) ??
                    'UNSPECIFIED',
                },
              ]}
              optionalFields={[
                { label: 'Keep Phase', render: item => item.keepPhase },
                {
                  label: 'Target PH',
                  render: item =>
                    item.targetPh != null && item.targetPh !== 0 ? item.targetPh : '',
                },
                {
                  label: 'Duration',
                  render: item =>
                    item.duration ? vpu(item.duration, ord.Time.TimeUnit) : '',
                },
                {
                  label: 'Automated',
                  render: item =>
                    item.isAutomated == null ? '' : reactionBoolean(item.isAutomated),
                },
                { label: 'Details', render: item => item.details },
                {
                  label: 'Aliqout Amount',
                  render: item => (item.amount ? formatAmount(item.amount) : ''),
                },
                {
                  label: 'Temperature',
                  render: item =>
                    item.temperature ? (
                      <RequiredOptionalFields
                        entity={item}
                        requiredFields={[
                          {
                            label: 'Setpoint',
                            render: ({ temperature }) =>
                              temperature?.setpoint
                                ? vpu(
                                    temperature.setpoint,
                                    ord.Temperature.TemperatureUnit,
                                  )
                                : 'UNSPECIFIED',
                          },
                        ]}
                        optionalFields={[
                          {
                            label: 'Control Type',
                            render: ({ temperature }) =>
                              enumName(
                                ord.TemperatureConditions.TemperatureControl
                                  .TemperatureControlType,
                                temperature?.control?.type,
                              ),
                          },
                          {
                            label: 'Control Details',
                            render: ({ temperature }) => temperature?.control?.details,
                          },
                        ]}
                      />
                    ) : (
                      ''
                    ),
                },
                {
                  label: 'Stirring',
                  render: item =>
                    item.stirring ? (
                      <RequiredOptionalFields
                        entity={item}
                        requiredFields={[]}
                        optionalFields={[
                          {
                            label: 'Method',
                            render: ({ stirring }) =>
                              enumName(
                                ord.StirringConditions.StirringMethodType,
                                stirring?.type,
                              ),
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
                            render: ({ stirring }) =>
                              stirring?.rate?.rpm != null && stirring.rate.rpm !== 0
                                ? stirring.rate.rpm
                                : '',
                          },
                        ]}
                      />
                    ) : (
                      ''
                    ),
                },
              ]}
            />
            {workup.input && (
              <>
                <InputsHeaderRow />
                <Accordion
                  variant="separated"
                  chevronPosition="left"
                  multiple={true}
                  defaultValue={['Input']}
                >
                  <InputComponentsListItem
                    input={workup.input}
                    name="Input"
                  />
                </Accordion>
              </>
            )}
          </Flex>
        ))}
      </Flex>
    </Flex>
  );
}
