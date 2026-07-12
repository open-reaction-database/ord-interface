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

// Ported from ord-app's KeyValueDisplay (common/components/display/KeyValueDisplay).
import { type ReactNode } from 'react';
import { Flex, Text, Tooltip } from '@mantine/core';
import clsx from 'clsx';
import typographyClasses from '../../styles/typography.module.scss';
import classes from './KeyValueDisplay.module.scss';

interface KeyValueDisplayProps {
  label: string;
  value: ReactNode;
  multiline?: boolean;
}

export function KeyValueDisplay({
  label,
  value,
  multiline,
}: Readonly<KeyValueDisplayProps>) {
  return (
    <Flex
      gap="xs"
      align="flex-start"
      className={clsx({ [classes.multilineWrapper]: !multiline })}
    >
      <Text className={clsx(typographyClasses.secondary2, classes.label)}>
        {label}:
      </Text>
      {multiline ? (
        <Text className={classes.multilineValue}>{value}</Text>
      ) : (
        <Tooltip label={value}>
          <Text className={classes.inlineValue}>{value}</Text>
        </Tooltip>
      )}
    </Flex>
  );
}

export default KeyValueDisplay;
