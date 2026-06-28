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

import React, { useState, useMemo } from 'react';
import { useLocation, useNavigate } from 'react-router-dom';
import SearchOptions from './SearchOptions';
import SearchResults from './SearchResults';
import LoadingSpinner from '../../components/LoadingSpinner';
import { useSearchTask } from '../../hooks/useSearchTask';
import './MainSearch.scss';

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
  const location = useLocation();
  const navigate = useNavigate();
  const [showOptions, setShowOptions] = useState(false);

  const hasSearchParams = useMemo(() => {
    const params = new URLSearchParams(location.search);
    return MEANINGFUL_QUERY_PARAMS.some(param => params.has(param));
  }, [location.search]);

  const { data, isFetching, error } = useSearchTask(location.search, hasSearchParams);
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

    navigate({ pathname: '/search', search: searchParams.toString() });
  };

  return (
    <div id="search-main">
      <div
        className={`search-options-container ${showOptions ? 'slide-out' : 'hidden'}`}
      >
        <div className="title">Filters & Options</div>
        <div className="options-holder">
          <SearchOptions onSearchOptions={updateSearchOptions} />
        </div>
        <div
          className="slide-out-tab"
          onClick={() => setShowOptions(!showOptions)}
        >
          <div className="line"></div>
          <div className="line"></div>
          <div className="line"></div>
        </div>
      </div>

      <div className="search-results">
        {!hasSearchParams && (
          <div className="no-results">
            <div className="title">
              Enter search criteria using the filters and options panel, then click
              Search.
            </div>
          </div>
        )}
        {hasSearchParams && error && (
          <div className="no-results">
            <div className="title">Search failed: {error.message}</div>
          </div>
        )}
        {hasSearchParams && !error && !loading && searchResults.length > 0 && (
          <SearchResults searchResults={searchResults} />
        )}
        {hasSearchParams && !error && !loading && !searchResults.length && (
          <div className="no-results">
            <div className="title">
              No results. Adjust the filters and options and search again.
            </div>
          </div>
        )}
        {hasSearchParams && !error && loading && (
          <div className="loading">
            <LoadingSpinner />
          </div>
        )}
      </div>
    </div>
  );
};

export default MainSearch;
