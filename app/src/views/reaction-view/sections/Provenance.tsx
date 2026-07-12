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

// Ported from ord-app's ReactionView/Provenance (read-only). Adds the
// city/DOI/patent/publication fields that ord-interface has always surfaced.
import { Flex, Title } from '@mantine/core';
import type { ord } from 'ord-schema-protobufjs';
import { RequiredOptionalFields } from '../../../components/display/RequiredOptionalFields';
import { formatDateToDisplay } from '../../../utils/formatting';
import { useReaction } from '../context';
import classes from './sections.module.scss';

export function Provenance() {
  const reaction = useReaction();
  const provenance: ord.IReactionProvenance = reaction.provenance ?? {};

  return (
    <Flex
      direction="column"
      gap="sm"
    >
      <Flex
        justify="space-between"
        align="center"
      >
        <Flex
          align="center"
          gap="sm"
        >
          <Title order={2}>Provenance</Title>
        </Flex>
      </Flex>
      <span className={classes.descriptionText}>
        Additional metadata about how this reaction was performed and originally
        reported
      </span>

      <Flex
        direction="column"
        gap="sm"
        className={classes.mainInformation}
      >
        <span className={classes.sectionLabel}>Experiment</span>
        <RequiredOptionalFields
          entity={provenance}
          requiredFields={[
            {
              label: 'Experimenter name',
              render: ({ experimenter }) => experimenter?.name,
            },
            { label: 'E-mail', render: ({ experimenter }) => experimenter?.email },
            { label: 'ORCID ID', render: ({ experimenter }) => experimenter?.orcid },
            {
              label: 'Username',
              render: ({ experimenter }) => experimenter?.username,
            },
          ]}
          optionalFields={[
            { label: 'City', render: p => p.city },
            {
              label: 'Experiment start',
              render: p => formatDateToDisplay(p.experimentStart),
            },
            { label: 'DOI', render: p => p.doi },
            { label: 'Patent', render: p => p.patent },
            {
              label: 'Publication URL',
              render: p =>
                p.publicationUrl ? (
                  <a
                    target="_blank"
                    href={p.publicationUrl}
                    rel="noreferrer"
                  >
                    {p.publicationUrl}
                  </a>
                ) : (
                  ''
                ),
            },
          ]}
        />

        <span className={classes.sectionLabel}>Record Creation</span>
        <RequiredOptionalFields
          entity={provenance.recordCreated ?? {}}
          requiredFields={[
            {
              label: 'Time',
              render: ({ time }) => (time ? formatDateToDisplay(time) : ''),
            },
            { label: 'E-mail', render: ({ person }) => person?.email },
            { label: 'ORCID ID', render: ({ person }) => person?.orcid },
            { label: 'Username', render: ({ person }) => person?.username },
            { label: 'Experimenter name', render: ({ person }) => person?.name },
          ]}
        />
      </Flex>

      <Flex
        direction="column"
        gap="sm"
      >
        {(provenance.recordModified ?? []).map((recordModification, index) => (
          <Flex
            direction="column"
            gap="xs"
            key={index}
          >
            <Flex
              align="center"
              gap="xs"
            >
              <Title order={3}>Record Modification {index + 1}</Title>
            </Flex>
            <RequiredOptionalFields
              entity={recordModification}
              requiredFields={[
                {
                  label: 'Time',
                  render: ({ time }) => (time ? formatDateToDisplay(time) : ''),
                },
                { label: 'Person Email', render: ({ person }) => person?.email },
                { label: 'Details', render: ({ details }) => details },
              ]}
            />
          </Flex>
        ))}
      </Flex>
    </Flex>
  );
}
