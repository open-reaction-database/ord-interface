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

// The shared component table row, ported from ord-app's ComponentsList /
// ComponentDisplayRowCustomActions (read-only: no action column contents).
/* eslint-disable react-refresh/only-export-components -- the grid classes are
   exported alongside the row component, mirroring ord-app's index. */
import { type ReactNode } from 'react';
import { Badge, Flex } from '@mantine/core';
import clsx from 'clsx';
import { ord } from 'ord-schema-protobufjs';
import { enumName } from '../../utils/enum';
import { KeyValueDisplay } from '../../components/display/KeyValueDisplay';
import CompoundPreview, { type ComponentLike } from './CompoundPreview';
import classes from './componentsList.module.scss';

export const componentsListClasses = classes;

interface ComponentDisplayRowProps {
  component: ComponentLike;
  renderDetails: (component: ComponentLike) => ReactNode;
  gridClassName?: string;
}

export function ComponentDisplayRow({
  component,
  renderDetails,
  gridClassName = clsx(classes.grid, classes.row),
}: Readonly<ComponentDisplayRowProps>) {
  const role = enumName(ord.ReactionRole.ReactionRoleType, component.reactionRole);

  return (
    <div className={gridClassName}>
      <Flex
        className={classes.identifiers}
        direction="column"
      >
        {(component.identifiers ?? []).map((identifier, index) => (
          <KeyValueDisplay
            key={index}
            label={
              enumName(
                ord.CompoundIdentifier.CompoundIdentifierType,
                identifier.type,
              ) ?? 'IDENTIFIER'
            }
            value={identifier.value}
          />
        ))}
      </Flex>
      <Flex
        align="center"
        justify="flex-start"
        className={clsx(classes.preview, classes.imagePreview)}
      >
        <CompoundPreview component={component} />
      </Flex>
      <Flex
        align="center"
        direction="column"
        gap={4}
        className={classes.role}
      >
        {role && role !== 'UNSPECIFIED' ? role : ''}
        {component.isLimiting === true && (
          <Badge
            size="xs"
            variant="light"
            w="fit-content"
          >
            Limiting reactant
          </Badge>
        )}
      </Flex>
      <Flex
        align="center"
        className={classes.details}
      >
        {renderDetails(component)}
      </Flex>
      <Flex
        className={classes.actions}
        align="center"
        justify="flex-end"
      ></Flex>
    </div>
  );
}
