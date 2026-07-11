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

import React from 'react';
import { useParams } from 'wouter';
import { Badge, Paper, SimpleGrid, Text, Title } from '@mantine/core';
import { useQuery } from '@tanstack/react-query';
import LoadingSpinner from '../../components/LoadingSpinner';
import CopyButton from '../../components/CopyButton';
import ChartView from './ChartView';
import SearchResults from './SearchResults';
import { useSearchTask } from '../../hooks/useSearchTask';
import { api } from '../../utils/api';
import type { Dataset } from '../../types/search';
import classes from './MainDatasetView.module.scss';

const MainDatasetView: React.FC = () => {
  const { datasetId } = useParams<{ datasetId: string }>();

  const datasetQuery = `?dataset_id=${datasetId}&limit=100`;
  const {
    data: searchData,
    isFetching,
    error,
  } = useSearchTask(datasetId ? datasetQuery : null, !!datasetId);
  const searchResults = searchData?.status === 'success' ? searchData.results : [];
  const loading = !!datasetId && (isFetching || searchData?.status === 'pending');

  const { data: datasetData, error: datasetError } = useQuery<Dataset>({
    queryKey: ['dataset-metadata', datasetId],
    enabled: !!datasetId,
    queryFn: () =>
      api
        .get<Dataset>('/dataset', { params: { dataset_id: datasetId } })
        .then(res => res.data),
  });

  return (
    <div className={classes.datasetMain}>
      {/* Dataset header, in the style of ord-app's DatasetHeader. */}
      <Paper
        radius="sm"
        p="lg"
        className={classes.header}
      >
        {datasetError ? (
          <Text c="secondary.1">
            Failed to load dataset metadata: {datasetError.message}
          </Text>
        ) : (
          <>
            <div className={classes.titleRow}>
              <Title order={1}>{datasetData?.name || datasetId}</Title>
              {datasetData && (
                <Badge
                  className={classes.countBadge}
                  variant="light"
                  color="gray"
                  radius="xl"
                >
                  {datasetData.num_reactions.toLocaleString()} reactions
                </Badge>
              )}
            </div>
            <div className={classes.datasetId}>
              {datasetId}
              <CopyButton textToCopy={datasetId ?? ''} />
            </div>
            {datasetData?.description && (
              <Text className={classes.description}>{datasetData.description}</Text>
            )}
          </>
        )}
      </Paper>

      {/* Composition charts. */}
      {datasetId && (
        <Paper
          radius="sm"
          p="lg"
          className={classes.charts}
        >
          <Title
            order={2}
            className={classes.chartsTitle}
          >
            Composition
          </Title>
          <SimpleGrid
            cols={{ base: 1, md: 2 }}
            spacing="lg"
          >
            <ChartView
              uniqueId="reactantsFrequency"
              title="Frequency of reactants"
              apiCall="input_stats"
              role="reactant"
              datasetId={datasetId}
            />
            <ChartView
              uniqueId="productsFrequency"
              title="Frequency of products"
              apiCall="product_stats"
              role="product"
              datasetId={datasetId}
            />
          </SimpleGrid>
        </Paper>
      )}

      {/* Reactions. */}
      <div className={classes.reactions}>
        {error ? (
          <Paper
            radius="sm"
            p="xl"
            className={classes.emptyState}
          >
            <Text c="secondary.1">Failed to load reactions: {error.message}</Text>
          </Paper>
        ) : loading ? (
          <LoadingSpinner />
        ) : searchResults.length > 0 ? (
          <SearchResults
            searchResults={searchResults}
            isOverflow={(datasetData?.num_reactions ?? 0) > searchResults.length}
          />
        ) : (
          <Paper
            radius="sm"
            p="xl"
            className={classes.emptyState}
          >
            <Text c="secondary.1">This dataset contains no reactions.</Text>
          </Paper>
        )}
      </div>
    </div>
  );
};

export default MainDatasetView;
