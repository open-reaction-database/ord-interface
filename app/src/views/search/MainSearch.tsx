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

import React, { useMemo } from 'react';
import { useLocation, useSearch } from 'wouter';
import { Paper, Text, Title } from '@mantine/core';
import SearchOptions from './SearchOptions';
import SearchResults from './SearchResults';
import LoadingSpinner from '../../components/LoadingSpinner';
import { useSearchTask } from '../../hooks/useSearchTask';
import classes from './MainSearch.module.scss';

interface SearchOptionsData {
  reagent: {
    reagents: Array<{ smileSmart: string; source: string; matchMode: string }>;
    useStereochemistry: boolean;
    similarityThreshold: number;
  };
  reaction: {
    reactionIds: string[];
    reactionSmarts: string[];
    min_yield: number;
    max_yield: number;
    min_conversion: number;
    max_conversion: number;
  };
  dataset: {
    datasetIds: string[];
    DOIs: string[];
  };
  general: {
    limit: number;
  };
}

const MEANINGFUL_QUERY_PARAMS = [
  'component',
  'dataset_id',
  'doi',
  'reaction_id',
  'reaction_smarts',
  'min_yield',
  'max_yield',
  'min_conversion',
  'max_conversion',
];

const MainSearch: React.FC = () => {
  const search = useSearch();
  const [, navigate] = useLocation();

  // The submit_query endpoint takes the page's query string verbatim.
  const queryString = search ? `?${search}` : '';

  const hasSearchParams = useMemo(() => {
    const params = new URLSearchParams(search);
    return MEANINGFUL_QUERY_PARAMS.some(param => params.has(param));
  }, [search]);

  const { data, isFetching, error } = useSearchTask(queryString, hasSearchParams);
  const searchResults = data?.status === 'success' ? data.results : [];
  const loading = hasSearchParams && (isFetching || data?.status === 'pending');

  const updateSearchOptions = (options: SearchOptionsData) => {
    const searchParams = new URLSearchParams();

    if (options.reagent.reagents.length) {
      options.reagent.reagents.forEach(reagent => {
        searchParams.append(
          'component',
          JSON.stringify({
            pattern: reagent.smileSmart,
            target: reagent.source,
            mode: reagent.matchMode.toLowerCase(),
          }),
        );
      });
      searchParams.set(
        'use_stereochemistry',
        options.reagent.useStereochemistry.toString(),
      );
      searchParams.set('similarity', options.reagent.similarityThreshold.toString());
    }

    options.dataset.datasetIds.forEach(id => searchParams.append('dataset_id', id));
    options.dataset.DOIs.forEach(doi => searchParams.append('doi', doi));
    options.reaction.reactionIds.forEach(id => searchParams.append('reaction_id', id));
    options.reaction.reactionSmarts.forEach(smarts =>
      searchParams.append('reaction_smarts', smarts),
    );

    // Yield/conversion: only include when narrower than the full 0-100 range.
    if (options.reaction.min_yield !== 0 || options.reaction.max_yield !== 100) {
      searchParams.set('min_yield', options.reaction.min_yield.toString());
      searchParams.set('max_yield', options.reaction.max_yield.toString());
    }
    if (
      options.reaction.min_conversion !== 0 ||
      options.reaction.max_conversion !== 100
    ) {
      searchParams.set('min_conversion', options.reaction.min_conversion.toString());
      searchParams.set('max_conversion', options.reaction.max_conversion.toString());
    }

    searchParams.set('limit', options.general.limit.toString() || '100');

    navigate(`/search?${searchParams.toString()}`);
  };

  return (
    <div className={classes.searchMain}>
      <aside className={classes.optionsPanel}>
        <Paper
          radius="sm"
          p="md"
          className={classes.optionsPaper}
        >
          <Title
            order={2}
            className={classes.optionsTitle}
          >
            Filters & options
          </Title>
          <SearchOptions onSearchOptions={updateSearchOptions} />
        </Paper>
      </aside>

      <div className={classes.results}>
        {!hasSearchParams && (
          <Paper
            radius="sm"
            p="xl"
            className={classes.emptyState}
          >
            <Text className={classes.emptyText}>
              Enter search criteria using the filters and options panel, then click
              Search.
            </Text>
          </Paper>
        )}
        {hasSearchParams && error && (
          <Paper
            radius="sm"
            p="xl"
            className={classes.emptyState}
          >
            <Text className={classes.emptyText}>Search failed: {error.message}</Text>
          </Paper>
        )}
        {hasSearchParams && !error && !loading && searchResults.length > 0 && (
          <SearchResults searchResults={searchResults} />
        )}
        {hasSearchParams && !error && !loading && !searchResults.length && (
          <Paper
            radius="sm"
            p="xl"
            className={classes.emptyState}
          >
            <Text className={classes.emptyText}>
              No results. Adjust the filters and options and search again.
            </Text>
          </Paper>
        )}
        {hasSearchParams && !error && loading && <LoadingSpinner />}
      </div>
    </div>
  );
};

export default MainSearch;
