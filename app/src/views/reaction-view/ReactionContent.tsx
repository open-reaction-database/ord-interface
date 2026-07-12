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

// Ported from ord-app's ReactionContent + ReactionTabs: the same nine section
// components rendered as tabs or stacked as a list.
import React, { Fragment, type FC } from 'react';
import { Stack, Tabs, Tooltip } from '@mantine/core';
import { Conditions } from './sections/Conditions';
import { Identifiers } from './sections/Identifiers';
import { Inputs } from './sections/Inputs';
import { Notes } from './sections/Notes';
import { Observations } from './sections/Observations';
import { Outcomes } from './sections/Outcomes';
import { Provenance } from './sections/Provenance';
import { Setup } from './sections/Setup';
import { Workups } from './sections/Workups';
import classes from './reactionTabs.module.scss';

interface ReactionTab {
  name: string;
  required?: true;
  Component: FC;
}

const tabs: Array<ReactionTab> = [
  { name: 'inputs', required: true, Component: Inputs },
  { name: 'outcomes', required: true, Component: Outcomes },
  { name: 'conditions', Component: Conditions },
  { name: 'identifiers', Component: Identifiers },
  { name: 'setup', Component: Setup },
  { name: 'notes', Component: Notes },
  { name: 'observations', Component: Observations },
  { name: 'workups', Component: Workups },
  { name: 'provenance', required: true, Component: Provenance },
];

function RequiredAsterisk() {
  return <span className={classes.required}>*</span>;
}

export function ReactionTabs() {
  return (
    <Tabs
      defaultValue={tabs[0].name}
      classNames={{ tab: classes.tabTitle, panel: classes.panel }}
    >
      <Tabs.List>
        {tabs.map(({ name, required }) => (
          <Fragment key={name}>
            {required ? (
              <Tooltip label="Mandatory section">
                <Tabs.Tab value={name}>
                  {name}
                  <RequiredAsterisk />
                </Tabs.Tab>
              </Tooltip>
            ) : (
              <Tabs.Tab value={name}>{name}</Tabs.Tab>
            )}
          </Fragment>
        ))}
      </Tabs.List>
      {tabs.map(({ name, Component }) => (
        <Tabs.Panel
          key={name}
          value={name}
        >
          <Component />
        </Tabs.Panel>
      ))}
    </Tabs>
  );
}

interface ReactionContentProps {
  viewMode: 'tabs' | 'list';
}

export const ReactionContent = React.memo(
  ({ viewMode }: Readonly<ReactionContentProps>) => {
    if (viewMode === 'tabs') {
      return <ReactionTabs />;
    }
    return (
      <Stack gap="xl">
        <Inputs />
        <Outcomes />
        <Conditions />
        <Identifiers />
        <Setup />
        <Notes />
        <Observations />
        <Workups />
        <Provenance />
      </Stack>
    );
  },
);

ReactionContent.displayName = 'ReactionContent';
