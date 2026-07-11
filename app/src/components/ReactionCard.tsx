/**
 * Copyright 2023 Open Reaction Database Project Authors
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

import React, { useCallback, useEffect, useState } from 'react';
import { Link, useLocation } from 'wouter';
import { Badge, Button, Checkbox, Loader, Paper, Text } from '@mantine/core';
import { ord } from 'ord-schema-protobufjs';
import clsx from 'clsx';
import CopyButton from './CopyButton';
import { api } from '../utils/api';
import { enumName } from '../utils/enum';
import { formatPercentage } from '../utils/outcomes';
import type { ReactionData, SearchResult } from '../types/search';
import classes from './ReactionCard.module.scss';

interface ReactionCardProps {
  reaction: SearchResult;
  isSelectable?: boolean;
  isSelected?: boolean;
  onSelectionChange?: (reactionId: string, isSelected: boolean) => void;
}

const YIELD_MEASUREMENT_TYPE = ord.ProductMeasurement.ProductMeasurementType.YIELD;

const ReactionCard: React.FC<ReactionCardProps> = ({
  reaction,
  isSelectable = true,
  isSelected = false,
  onSelectionChange,
}) => {
  const [, navigate] = useLocation();
  const [reactionTable, setReactionTable] = useState<string | null>(null);
  const [tableError, setTableError] = useState(false);

  const getReactionTable = useCallback(async () => {
    try {
      const response = await api.get<string>('/reaction_summary', {
        params: { reaction_id: reaction.reaction_id },
      });
      setReactionTable(response.data);
    } catch (error) {
      // Never pipe an HTML error body into dangerouslySetInnerHTML.
      console.error('Error fetching reaction table:', error);
      setTableError(true);
    }
  }, [reaction.reaction_id]);

  const getYield = (measurements: ord.IProductMeasurement[]): string => {
    const yieldObj = measurements.find(m => m.type === YIELD_MEASUREMENT_TYPE);
    return yieldObj?.percentage ? formatPercentage(yieldObj.percentage) : 'Not listed';
  };

  const getConversion = (data: ReactionData | undefined): string => {
    const conversion = data?.outcomes?.[0]?.conversion;
    if (!conversion) return 'Not listed';
    return formatPercentage(conversion);
  };

  const conditionsAndDuration = (data: ReactionData | undefined): string[] => {
    const details: string[] = [];
    if (!data) return details;

    const temp = data.conditions?.temperature?.setpoint;
    if (temp) {
      const units = enumName(ord.Temperature.TemperatureUnit, temp.units);
      details.push(`at ${temp.value ?? 0}${units ? ` ${units.toLowerCase()}` : '°C'}`);
    }

    const pressure = data.conditions?.pressure?.setpoint;
    if (pressure) {
      const units = enumName(ord.Pressure.PressureUnit, pressure.units);
      details.push(
        `under ${pressure.value ?? 0}${units ? ` ${units.toLowerCase()}` : ' atm'}`,
      );
    }

    const reactionTime = data.outcomes?.[0]?.reactionTime;
    if (reactionTime?.value) {
      const units = enumName(ord.Time.TimeUnit, reactionTime.units);
      details.push(
        `for ${reactionTime.value}${units ? ` ${units.toLowerCase()}` : 's'}`,
      );
    }

    return details;
  };

  const productIdentifier = (identifier: ord.ICompoundIdentifier): string => {
    const type = enumName(
      ord.CompoundIdentifier.CompoundIdentifierType,
      identifier.type,
    );
    return `${type ?? ''}: ${identifier.value ?? ''}`;
  };

  useEffect(() => {
    getReactionTable();
  }, [getReactionTable]);

  const reactionData = reaction.data;
  const firstOutcome = reactionData?.outcomes?.[0];
  const firstProduct = firstOutcome?.products?.[0];
  const firstProductIdentifier = firstProduct?.identifiers?.[0];
  const provenance = reactionData?.provenance;

  return (
    <Paper
      radius="sm"
      p="lg"
      className={clsx(classes.card, { [classes.selected]: isSelected })}
    >
      <div className={classes.topRow}>
        <div className={classes.topLeft}>
          {isSelectable && (
            <Checkbox
              size="sm"
              checked={isSelected}
              onChange={event =>
                onSelectionChange?.(reaction.reaction_id, event.currentTarget.checked)
              }
              aria-label="Select reaction"
            />
          )}
          <Text
            className={classes.reactionId}
            title={reaction.reaction_id}
          >
            {reaction.reaction_id}
          </Text>
          {provenance?.isMined && (
            <Badge
              variant="light"
              color="orange"
              radius="xl"
              className={classes.minedBadge}
            >
              Mined
            </Badge>
          )}
        </div>
        <Button
          variant="default"
          size="compact-sm"
          radius="sm"
          onClick={() => navigate(`/id/${reaction.reaction_id}`)}
        >
          View details
        </Button>
      </div>

      <div className={classes.preview}>
        {reactionTable ? (
          <div dangerouslySetInnerHTML={{ __html: reactionTable }} />
        ) : tableError ? (
          <Text className={classes.previewError}>Reaction preview unavailable.</Text>
        ) : (
          <Loader size="sm" />
        )}
      </div>

      <div className={classes.descriptors}>
        <div className={classes.descriptorColumn}>
          <div className={classes.descriptor}>
            <span className={classes.label}>Yield</span>
            <span>{getYield(firstProduct?.measurements ?? [])}</span>
          </div>
          <div className={classes.descriptor}>
            <span className={classes.label}>Conversion</span>
            <span>{getConversion(reactionData)}</span>
          </div>
          <div className={classes.descriptor}>
            <span className={classes.label}>Conditions</span>
            <span>
              {conditionsAndDuration(reactionData).join('; ') || 'Not listed'}
            </span>
          </div>
          {firstProductIdentifier && (
            <div className={classes.descriptor}>
              <span className={classes.label}>Product</span>
              <span className={classes.productValue}>
                {productIdentifier(firstProductIdentifier)}
                <CopyButton textToCopy={firstProductIdentifier.value || ''} />
              </span>
            </div>
          )}
        </div>

        <div className={classes.descriptorColumn}>
          <div className={classes.descriptor}>
            <span className={classes.label}>Uploaded by</span>
            <span>
              {provenance?.recordCreated?.person?.name || 'Unknown'},{' '}
              {provenance?.recordCreated?.person?.organization || 'Unknown'}
            </span>
          </div>
          <div className={classes.descriptor}>
            <span className={classes.label}>Uploaded on</span>
            <span>
              {provenance?.recordCreated?.time?.value
                ? new Date(provenance.recordCreated.time.value).toLocaleDateString()
                : 'Unknown'}
            </span>
          </div>
          <div className={classes.descriptor}>
            <span className={classes.label}>DOI</span>
            <span>{provenance?.doi || 'Not available'}</span>
          </div>
          {provenance?.publicationUrl && (
            <div className={classes.descriptor}>
              <span className={classes.label}>Publication</span>
              <a
                href={provenance.publicationUrl}
                target="_blank"
                rel="noopener noreferrer"
              >
                {provenance.publicationUrl}
              </a>
            </div>
          )}
          {reaction.dataset_id && (
            <div className={classes.descriptor}>
              <span className={classes.label}>Dataset</span>
              <Link href={`/dataset/${reaction.dataset_id}`}>
                {reaction.dataset_id}
              </Link>
            </div>
          )}
        </div>
      </div>
    </Paper>
  );
};

export default ReactionCard;
