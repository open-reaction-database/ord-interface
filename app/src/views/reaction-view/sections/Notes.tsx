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

// Ported from ord-app's ReactionView/Notes (read-only): a 2-column grid of
// whichever note fields are set.
import { Fragment, useMemo } from 'react';
import { Flex, Title } from '@mantine/core';
import type { ord } from 'ord-schema-protobufjs';
import typographyClasses from '../../../styles/typography.module.scss';
import { reactionBoolean } from '../../../utils/formatting';
import { useReaction } from '../context';
import classes from './sections.module.scss';

const NOTES_FIELDS: Array<[keyof ord.IReactionNotes, string]> = [
  ['procedureDetails', 'Procedure details'],
  ['safetyNotes', 'Safety notes'],
  ['isHeterogeneous', 'Is heterogeneous'],
  ['formsPrecipitate', 'Forms precipitate'],
  ['isExothermic', 'Is exothermic'],
  ['offgasses', 'Offgasses'],
  ['isSensitiveToOxygen', 'Oxygen sensitive'],
  ['isSensitiveToMoisture', 'Moisture sensitive'],
  ['isSensitiveToLight', 'Light sensitive'],
];

export function Notes() {
  const reaction = useReaction();

  const fields = useMemo((): Array<[string, string]> => {
    const notes: ord.IReactionNotes = reaction.notes ?? {};
    return NOTES_FIELDS.map(([key, label]): [string, string] => {
      const value = notes[key];
      if (typeof value === 'string') return [label, value];
      if (typeof value === 'boolean') return [label, reactionBoolean(value)];
      return [label, ''];
    }).filter(([, value]) => value && value !== 'Unspecified');
  }, [reaction.notes]);

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
          <Title order={2}>Notes</Title>
        </Flex>
      </Flex>
      <div className={classes.notesGrid}>
        {fields.map(([label, value]) => (
          <Fragment key={label}>
            <span className={typographyClasses.secondary2}>{label}</span>
            <p>{value}</p>
          </Fragment>
        ))}
      </div>
    </Flex>
  );
}
