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

// Ported from ord-app's ReactionView/Observation (read-only).
import React from 'react';
import { Flex, Title } from '@mantine/core';
import { ord } from 'ord-schema-protobufjs';
import { Counter } from '../../../components/display/Counter';
import { RequiredOptionalFields } from '../../../components/display/RequiredOptionalFields';
import { vpu } from '../../../utils/formatting';
import { useReaction } from '../context';

// Renders the displayable slice of a Data message (ord-app's AppDataDisplay).
function dataDisplay(data: ord.IData | null | undefined): React.ReactNode {
  if (!data) return '';
  if (data.url) {
    return (
      <a
        target="_blank"
        href={data.url}
        rel="noreferrer"
      >
        {data.url}
      </a>
    );
  }
  if (data.stringValue) return data.stringValue;
  if (data.floatValue != null) return String(data.floatValue);
  if (data.integerValue != null) return String(data.integerValue);
  if (data.bytesValue?.length) {
    return `(binary data, ${data.bytesValue.length.toLocaleString()} bytes)`;
  }
  return '';
}

export function Observations() {
  const reaction = useReaction();
  const observations = reaction.observations ?? [];

  return (
    <Flex direction="column">
      <Flex justify="space-between">
        <Flex
          align="center"
          gap="sm"
        >
          <Title order={2}>Observations</Title>
          <Counter amount={observations.length} />
        </Flex>
      </Flex>
      <span>
        Observations are time-stamped comments, images, etc. that are recorded during
        the reaction
      </span>
      <Flex
        direction="column"
        gap="sm"
      >
        {observations.map((observation, index) => (
          <Flex
            direction="column"
            gap="xs"
            key={index}
          >
            <Flex
              align="center"
              gap="xs"
            >
              <Title order={3}>Observation {index + 1}</Title>
            </Flex>
            <RequiredOptionalFields
              entity={observation}
              requiredFields={[
                {
                  label: 'Time',
                  render: ({ time }) =>
                    time?.value ? vpu(time, ord.Time.TimeUnit) : '',
                },
              ]}
              optionalFields={[
                { label: 'Comment', render: ({ comment }) => comment },
                { label: 'Data', render: ({ image }) => dataDisplay(image) },
                {
                  label: 'Description',
                  render: ({ image }) => image?.description,
                },
              ]}
            />
          </Flex>
        ))}
      </Flex>
    </Flex>
  );
}
