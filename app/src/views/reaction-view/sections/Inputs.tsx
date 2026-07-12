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

// Ported from ord-app's ReactionView/Inputs (read-only): header + a separated
// accordion of inputs, each an identifiers/preview/role/amount grid.
import { useMemo } from 'react';
import { Accordion, Flex, Text, Title } from '@mantine/core';
import { IconInbox } from '@tabler/icons-react';
import clsx from 'clsx';
import { ord } from 'ord-schema-protobufjs';
import { Counter } from '../../../components/display/Counter';
import typographyClasses from '../../../styles/typography.module.scss';
import { formatAmount } from '../../../utils/formatting';
import { sortedInputEntries } from '../../../utils/proto';
import { useReaction } from '../context';
import { ComponentDisplayRow, componentsListClasses } from '../ComponentsList';
import type { ComponentLike } from '../CompoundPreview';
import classes from './inputs.module.scss';

const renderDetails = (component: ComponentLike) => formatAmount(component.amount);

interface InputItemProps {
  name: string;
  input: ord.IReactionInput;
}

// Exported for reuse by the Workups section (workup.input rendering).
export function InputComponentsListItem({ name, input }: Readonly<InputItemProps>) {
  const components = input.components ?? [];

  return (
    <Accordion.Item value={name}>
      <Accordion.Control className={classes.control}>{name}</Accordion.Control>
      <Accordion.Panel>
        {components.length > 0 ? (
          components.map((component, index) => (
            <ComponentDisplayRow
              key={index}
              component={component as ComponentLike}
              renderDetails={renderDetails}
              gridClassName={clsx(classes.grid, classes.row)}
            />
          ))
        ) : (
          <Flex
            direction="column"
            align="center"
            gap="8"
          >
            <Text className={typographyClasses.secondary1}>
              There are no Components yet
            </Text>
          </Flex>
        )}
      </Accordion.Panel>
    </Accordion.Item>
  );
}

export function InputsHeaderRow() {
  const headers = [
    { label: 'Input', className: classes.input },
    { label: 'Identifiers', className: componentsListClasses.identifiers },
    { label: 'Preview', className: componentsListClasses.preview },
    { label: 'Role', className: componentsListClasses.role },
    { label: 'Amount', className: componentsListClasses.details },
  ];
  return (
    <div className={classes.grid}>
      {headers.map(({ label, className }) => (
        <Text
          size="md"
          className={clsx(classes.text, className)}
          key={label}
        >
          {label}
        </Text>
      ))}
    </div>
  );
}

export function Inputs() {
  const reaction = useReaction();
  const inputs = useMemo(() => sortedInputEntries(reaction.inputs), [reaction]);
  const names = inputs.map(([name]) => name);

  return (
    <Flex direction="column">
      <Flex justify="space-between">
        <Flex
          align="center"
          gap="sm"
        >
          <Title order={2}>Inputs</Title>
          <Counter amount={inputs.length} />
        </Flex>
      </Flex>
      <span>Reaction inputs include every chemical added to the reaction vessel</span>
      {inputs.length > 0 ? (
        <>
          <InputsHeaderRow />
          <Accordion
            variant="separated"
            chevronPosition="left"
            multiple={true}
            defaultValue={names}
          >
            {inputs.map(([name, input]) => (
              <InputComponentsListItem
                key={name}
                name={name}
                input={input}
              />
            ))}
          </Accordion>
        </>
      ) : (
        <Flex
          direction="column"
          align="center"
          gap="sm"
        >
          <IconInbox
            className={classes.icon}
            stroke={1}
          />
          <span className={typographyClasses.secondary1}>There are no Inputs yet</span>
        </Flex>
      )}
    </Flex>
  );
}
