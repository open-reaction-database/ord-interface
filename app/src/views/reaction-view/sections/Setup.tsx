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

// Ported from ord-app's ReactionView/Setup (read-only).
import { Flex, Title } from '@mantine/core';
import { ord } from 'ord-schema-protobufjs';
import { RequiredOptionalFields } from '../../../components/display/RequiredOptionalFields';
import { enumName } from '../../../utils/enum';
import { reactionBoolean, vpu } from '../../../utils/formatting';
import { useReaction } from '../context';
import classes from './sections.module.scss';

export function Setup() {
  const reaction = useReaction();
  const setup: ord.IReactionSetup = reaction.setup ?? {};
  const vessel = setup.vessel ?? {};

  return (
    <Flex direction="column">
      <Flex justify="space-between">
        <Flex
          align="center"
          gap="sm"
        >
          <Title order={2}>Setup</Title>
        </Flex>
      </Flex>
      <Flex direction="column">
        <Flex className={classes.setupContainer}>
          <RequiredOptionalFields
            entity={vessel}
            requiredFields={[
              {
                label: 'Vessel',
                render: v => enumName(ord.Vessel.VesselType, v.type) ?? 'UNSPECIFIED',
              },
              { label: 'Details', render: v => v.details },
            ]}
          />
        </Flex>
        <RequiredOptionalFields
          entity={vessel}
          requiredFields={[
            {
              label: 'Material',
              render: v =>
                enumName(ord.VesselMaterial.VesselMaterialType, v.material?.type) ??
                'UNSPECIFIED',
            },
            { label: 'Details', render: v => v.material?.details },
          ]}
          optionalFields={[
            {
              label: 'Volume',
              render: v => (v.volume ? vpu(v.volume, ord.Volume.VolumeUnit) : ''),
            },
          ]}
        />
        <RequiredOptionalFields
          entity={setup}
          requiredFields={[]}
          optionalFields={[
            {
              label: 'Automated',
              render: s =>
                s.isAutomated == null ? '' : reactionBoolean(s.isAutomated),
            },
            { label: 'Automation platform', render: s => s.automationPlatform },
            {
              label: 'Environment',
              render: s =>
                enumName(
                  ord.ReactionSetup.ReactionEnvironment.ReactionEnvironmentType,
                  s.environment?.type,
                ),
            },
            { label: 'Environment details', render: s => s.environment?.details },
          ]}
        />
        {(vessel.preparations ?? []).map((preparation, index) => (
          <Flex
            direction="column"
            gap="xs"
            key={`prep-${index}`}
          >
            <Flex
              align="center"
              gap="xs"
            >
              <Title order={3}>Vessel Preparation {index + 1}</Title>
            </Flex>
            <RequiredOptionalFields
              entity={preparation}
              requiredFields={[
                {
                  label: 'Type',
                  render: ({ type }) =>
                    enumName(ord.VesselPreparation.VesselPreparationType, type) ??
                    'UNSPECIFIED',
                },
                { label: 'Details', render: ({ details }) => details },
              ]}
            />
          </Flex>
        ))}
        {(vessel.attachments ?? []).map((attachment, index) => (
          <Flex
            direction="column"
            gap="xs"
            key={`attach-${index}`}
          >
            <Flex
              align="center"
              gap="xs"
            >
              <Title order={3}>Vessel Attachments {index + 1}</Title>
            </Flex>
            <RequiredOptionalFields
              entity={attachment}
              requiredFields={[
                {
                  label: 'Type',
                  render: ({ type }) =>
                    enumName(ord.VesselAttachment.VesselAttachmentType, type) ??
                    'UNSPECIFIED',
                },
                { label: 'Details', render: ({ details }) => details },
              ]}
            />
          </Flex>
        ))}
      </Flex>
    </Flex>
  );
}
