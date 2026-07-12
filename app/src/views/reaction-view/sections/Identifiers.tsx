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

// Ported from ord-app's ReactionView/Identifiers (read-only).
import { Flex, Title } from '@mantine/core';
import { ord } from 'ord-schema-protobufjs';
import { Counter } from '../../../components/display/Counter';
import { KeyValueDisplay } from '../../../components/display/KeyValueDisplay';
import { enumName } from '../../../utils/enum';
import { useReaction } from '../context';
import classes from './sections.module.scss';

export function Identifiers() {
  const reaction = useReaction();
  const identifiers = reaction.identifiers ?? [];

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
          <Title order={2}>Identifiers</Title>
          <Counter amount={identifiers.length} />
        </Flex>
      </Flex>
      <span className={classes.descriptionText}>
        Reaction identifiers define descriptions of the overall reaction
      </span>
      <Flex
        direction="column"
        gap="sm"
      >
        {identifiers.map((identifier, index) => (
          <div key={index}>
            <Flex align="center">
              <span className={classes.sectionLabel}>Identifier {index + 1}</span>
            </Flex>
            <KeyValueDisplay
              label={
                enumName(
                  ord.ReactionIdentifier.ReactionIdentifierType,
                  identifier.type,
                ) ?? 'TYPE'
              }
              value={identifier.value}
              multiline
            />
            <KeyValueDisplay
              label="Details"
              value={identifier.details}
              multiline
            />
          </div>
        ))}
      </Flex>
    </Flex>
  );
}
