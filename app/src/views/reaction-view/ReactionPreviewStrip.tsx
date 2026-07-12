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

// The reaction "scheme" strip in the page header, ported from ord-app's
// ReactionPreview / ReactionInputPreview / ReactionOutcomePreview /
// ComponentMetadata, adapted to raw ord.I* messages.
import React, { Fragment, useMemo } from 'react';
import { Badge, Flex, Text, Tooltip } from '@mantine/core';
import { ord } from 'ord-schema-protobufjs';
import { enumName } from '../../utils/enum';
import { formatAmount, renderValuePrecisionUnit, vpu } from '../../utils/formatting';
import { sortedInputEntries } from '../../utils/proto';
import CompoundPreview, { type ComponentLike } from './CompoundPreview';
import classes from './reactionPreview.module.scss';

const NAME_TYPE = ord.CompoundIdentifier.CompoundIdentifierType.NAME;
const YIELD_TYPE = ord.ProductMeasurement.ProductMeasurementType.YIELD;

function getProductYieldPercent(product: ord.IProductCompound): number | undefined {
  const yieldValue = product.measurements?.find(
    measurement => measurement.type === YIELD_TYPE,
  )?.percentage?.value;
  return yieldValue ?? undefined;
}

const ComponentMetadata: React.FC<{ component: ComponentLike }> = ({ component }) => {
  const name = (component.identifiers ?? []).find(
    identifier => identifier.type === NAME_TYPE,
  );
  const productYield =
    'measurements' in component && component.measurements
      ? getProductYieldPercent(component)
      : undefined;
  const amount = component.amount ? formatAmount(component.amount) : '';
  const role = enumName(ord.ReactionRole.ReactionRoleType, component.reactionRole);

  return (
    <Flex
      direction="column"
      justify="flex-end"
      className={classes.componentsMetadata}
    >
      {name?.value && (
        <Tooltip label={name.value}>
          <Text
            size="xs"
            className={classes.name}
          >
            {name.value}
          </Text>
        </Tooltip>
      )}
      {amount && <Text size="xs">{amount}</Text>}
      {productYield != null && <Text size="xs">{productYield}% yield</Text>}
      {role && role !== 'UNSPECIFIED' && <Text size="xs">{role}</Text>}
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
  );
};

const ComponentCell: React.FC<{ component: ComponentLike }> = ({ component }) => (
  <div className={classes.component}>
    <div className={classes.molecule}>
      <CompoundPreview component={component} />
    </div>
    {'isDesiredProduct' in component && component.isDesiredProduct === true && (
      <Badge
        classNames={{ root: classes.desiredProductBadge, label: classes.badgeLabel }}
      >
        ✨ Desired
      </Badge>
    )}
    <ComponentMetadata component={component} />
  </div>
);

const PreviewCard: React.FC<{
  label: string;
  components: ComponentLike[];
}> = ({ label, components }) => (
  <div className={classes.inputCard}>
    <Badge
      variant="light"
      color="primary"
      size="lg"
    >
      {label}
    </Badge>
    <Flex
      gap="sm"
      align="center"
      className={classes.componentList}
    >
      {components.length === 0 ? (
        <CompoundPreview component={null} />
      ) : (
        components.map((component, index) => (
          <ComponentCell
            key={index}
            component={component}
          />
        ))
      )}
    </Flex>
  </div>
);

const Arrow: React.FC = () => (
  <svg
    className={classes.arrow}
    viewBox="0 0 120 24"
    fill="none"
    xmlns="http://www.w3.org/2000/svg"
  >
    <path
      d="M2 12h108m0 0-9-6m9 6-9 6"
      stroke="currentColor"
      strokeWidth="2"
      strokeLinecap="round"
      strokeLinejoin="round"
    />
  </svg>
);

interface ReactionPreviewStripProps {
  reaction: ord.Reaction;
}

export const ReactionPreviewStrip: React.FC<ReactionPreviewStripProps> = ({
  reaction,
}) => {
  const inputs = useMemo(() => sortedInputEntries(reaction.inputs), [reaction]);
  const outcomes = reaction.outcomes ?? [];

  const outcomeLabel = (outcome: ord.IReactionOutcome): string => {
    const time = outcome.reactionTime?.value
      ? vpu(outcome.reactionTime, ord.Time.TimeUnit)
      : '';
    const conversion =
      outcome.conversion?.value != null
        ? `${renderValuePrecisionUnit({
            value: outcome.conversion.value,
            units: '%',
          })} conversion`
        : '';
    const parts = [time, conversion].filter(Boolean);
    return parts.length > 0 ? ` (${parts.join(', ')})` : '';
  };

  if (inputs.length === 0 && outcomes.length === 0) {
    return null;
  }

  return (
    <div className={classes.wrapper}>
      {inputs.map(([name, input], index) => (
        <Fragment key={name}>
          {index > 0 && <span className={classes.plus}>+</span>}
          <PreviewCard
            label={name}
            components={(input.components ?? []) as ComponentLike[]}
          />
        </Fragment>
      ))}
      <Flex
        direction="column"
        className={classes.arrowWrapper}
        justify="center"
      >
        <Arrow />
      </Flex>
      {outcomes.map((outcome, index) => (
        <PreviewCard
          key={index}
          label={`Outcome${outcomeLabel(outcome)}`}
          components={(outcome.products ?? []) as ComponentLike[]}
        />
      ))}
    </div>
  );
};

export default ReactionPreviewStrip;
