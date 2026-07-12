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

// Ported from ord-app's ReactionView/Outcomes (+ OutcomeListItem,
// OutcomeListItemHeader, MeasurementsPreview), read-only over ord.I* messages.
import React from 'react';
import { Accordion, Flex, Text, Title, Tooltip } from '@mantine/core';
import { IconClock, IconFlask, IconInbox } from '@tabler/icons-react';
import clsx from 'clsx';
import { ord } from 'ord-schema-protobufjs';
import { Counter } from '../../../components/display/Counter';
import { KeyValueDisplay } from '../../../components/display/KeyValueDisplay';
import typographyClasses from '../../../styles/typography.module.scss';
import { enumName } from '../../../utils/enum';
import { formatAmount, renderValuePrecisionUnit, vpu } from '../../../utils/formatting';
import { useReaction } from '../context';
import { ComponentDisplayRow, componentsListClasses } from '../ComponentsList';
import type { ComponentLike } from '../CompoundPreview';
import classes from './outcomes.module.scss';

const measurementValue = (measurement: ord.IProductMeasurement): string => {
  if (measurement.percentage) {
    return renderValuePrecisionUnit({
      value: measurement.percentage.value,
      precision: measurement.percentage.precision,
      units: '%',
    });
  }
  if (measurement.amount) return formatAmount(measurement.amount);
  if (measurement.floatValue) {
    return renderValuePrecisionUnit({
      value: measurement.floatValue.value,
      precision: measurement.floatValue.precision,
    });
  }
  if (measurement.stringValue) return measurement.stringValue;
  return '';
};

const MeasurementsPreview: React.FC<{ product: ord.IProductCompound }> = ({
  product,
}) => (
  <Flex
    direction="column"
    gap="xs"
  >
    {(product.measurements ?? []).map((measurement, index) => {
      const type = enumName(
        ord.ProductMeasurement.ProductMeasurementType,
        measurement.type,
      );
      const value = measurementValue(measurement);
      return (
        <Flex
          key={index}
          gap="sm"
          wrap="nowrap"
          className={classes.measurementWrapper}
        >
          {type && <Text className={classes.measurementKeyType}>{type}</Text>}
          {measurement.analysisKey && (
            <Text className={classes.measurementKeyType}>
              {measurement.analysisKey}
            </Text>
          )}
          {value && (
            <Tooltip label={value}>
              <Text className={classes.measurementValue}>{value}</Text>
            </Tooltip>
          )}
        </Flex>
      );
    })}
  </Flex>
);

const renderDetails = (component: ComponentLike) => (
  <MeasurementsPreview product={component} />
);

const productHeaders = [
  { label: 'Identifiers', className: componentsListClasses.identifiers },
  { label: 'Preview', className: componentsListClasses.preview },
  { label: 'Role', className: componentsListClasses.role },
  { label: 'Measurements', className: componentsListClasses.details },
];

const OutcomeListItem: React.FC<{
  outcome: ord.IReactionOutcome;
  index: number;
}> = ({ outcome, index }) => {
  const analyses = Object.entries(outcome.analyses ?? {});
  const products = outcome.products ?? [];
  const time = outcome.reactionTime?.value
    ? vpu(outcome.reactionTime, ord.Time.TimeUnit)
    : '';
  const conversion =
    outcome.conversion?.value != null
      ? renderValuePrecisionUnit({
          value: outcome.conversion.value,
          precision: outcome.conversion.precision,
          units: '%',
        })
      : '';

  return (
    <Accordion.Item value={`outcome-${index}`}>
      <Accordion.Control classNames={{ label: classes.label }}>
        <Flex
          align="center"
          gap="xs"
        >
          <Title order={3}>Outcome</Title>
          <span>·</span>
          {products.length}
        </Flex>
        {time && (
          <Flex
            align="center"
            gap="xs"
          >
            <IconClock className={classes.shortInfoIcon} />
            <Text className={clsx(classes.shortInfoText, typographyClasses.secondary2)}>
              Time:{' '}
            </Text>
            <Text className={classes.shortInfoText}>{time}</Text>
          </Flex>
        )}
        {conversion && (
          <Flex
            align="center"
            gap="xs"
          >
            <IconFlask className={classes.shortInfoIcon} />
            <Text className={clsx(classes.shortInfoText, typographyClasses.secondary2)}>
              Limiting reactant conversion:{' '}
            </Text>
            <Text className={classes.shortInfoText}>{conversion}</Text>
          </Flex>
        )}
      </Accordion.Control>
      <Accordion.Panel>
        {analyses.length > 0 && (
          <div className={classes.analysesList}>
            {analyses.map(([name, analysis]) => (
              <Flex
                key={name}
                className={classes.analysesItem}
                direction="column"
              >
                <KeyValueDisplay
                  label="Type"
                  value={enumName(ord.Analysis.AnalysisType, analysis.type)}
                />
                <KeyValueDisplay
                  label="Name"
                  value={name}
                />
                <KeyValueDisplay
                  label="Details"
                  value={analysis.details}
                />
              </Flex>
            ))}
          </div>
        )}
        <div className={componentsListClasses.grid}>
          {productHeaders.map(({ label, className }) => (
            <Text
              className={clsx(classes.text, className)}
              key={label}
            >
              {label}
            </Text>
          ))}
        </div>
        {products.map((product, productIndex) => (
          <ComponentDisplayRow
            key={productIndex}
            component={product as ComponentLike}
            renderDetails={renderDetails}
          />
        ))}
      </Accordion.Panel>
    </Accordion.Item>
  );
};

export function Outcomes() {
  const reaction = useReaction();
  const outcomes = reaction.outcomes ?? [];
  const ids = outcomes.map((_, index) => `outcome-${index}`);

  return (
    <Flex direction="column">
      <Flex justify="space-between">
        <Flex
          align="center"
          gap="sm"
        >
          <Title order={2}>Outcomes</Title>
          <Counter amount={outcomes.length} />
        </Flex>
      </Flex>
      <span>
        Outcomes record timestamped analyses and, optionally, product characterization
      </span>
      {outcomes.length > 0 ? (
        <Accordion
          variant="separated"
          chevronPosition="left"
          multiple={true}
          className={classes.itemsList}
          defaultValue={ids}
        >
          {outcomes.map((outcome, index) => (
            <OutcomeListItem
              key={index}
              outcome={outcome}
              index={index}
            />
          ))}
        </Accordion>
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
          <span className={typographyClasses.secondary1}>
            There are no Outcomes yet
          </span>
        </Flex>
      )}
    </Flex>
  );
}
