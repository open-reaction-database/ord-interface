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

import React, { useState } from 'react';
import { useParams } from 'react-router-dom';
import { useQuery } from '@tanstack/react-query';
import LoadingSpinner from '../../components/LoadingSpinner';
import ChartView from './ChartView';
import SearchResults from './SearchResults';
import { useSearchTask } from '../../hooks/useSearchTask';
import type { Dataset } from '../../types/search';
import './MainDatasetView.scss';

const MainDatasetView: React.FC = () => {
  const { datasetId } = useParams<{ datasetId: string }>();
  const [isCollapsed, setIsCollapsed] = useState(true);

  const datasetQuery = `?dataset_id=${datasetId}&limit=100`;
  const { data: searchData, isFetching, error } = useSearchTask(datasetId ? datasetQuery : null, !!datasetId);
  const searchResults = searchData?.status === 'success' ? searchData.results : [];
  const loading = !!datasetId && (isFetching || searchData?.status === 'pending');

  const { data: datasetData } = useQuery<Dataset | null>({
    queryKey: ['dataset-metadata', datasetId],
    enabled: !!datasetId,
    queryFn: async () => {
      const res = await fetch('/api/datasets');
      const datasets = (await res.json()) as Dataset[];
      return datasets.find(d => d.dataset_id === datasetId) ?? null;
    },
  });

  return (
    <div id="dataset-main">
      <h1>Dataset View</h1>
      <div
        id="charts"
        style={{
          display: isCollapsed ? 'grid' : 'flex',
          flexDirection: isCollapsed ? undefined : 'column',
        }}
      >
        <div
          id="chartsection"
          className="charts-container"
          style={{ width: isCollapsed ? '80%' : '100%' }}
        >
          <div className="charts-header">
            <div id="expand">
              <button onClick={() => setIsCollapsed(c => !c)}>
                <i
                  className="material-icons"
                  title={isCollapsed ? 'Expand' : 'Collapse'}
                >
                  {isCollapsed ? 'keyboard_double_arrow_right' : 'keyboard_double_arrow_left'}
                </i>
              </button>
            </div>
          </div>
          <div
            id="chartsectioncharts"
            className={`charts-content ${isCollapsed ? '' : 'expanded'}`}
          >
            <ChartView
              uniqueId="reactantsFrequency"
              title="Frequency of Reactants"
              apiCall="input_stats"
              role="reactant"
              isCollapsed={isCollapsed}
            />
            <ChartView
              uniqueId="productsFrequency"
              title="Frequency of Products"
              apiCall="product_stats"
              role="product"
              isCollapsed={isCollapsed}
            />
          </div>
        </div>

        <div id="datasection">
          <div className="h4">Dataset Metadata</div>
          <table>
            <tbody>
              <tr>
                <td>Dataset ID:</td>
                <td>{datasetData?.dataset_id ?? '(no id)'}</td>
              </tr>
              <tr>
                <td>Dataset Name:</td>
                <td>{datasetData?.name ?? '(no name)'}</td>
              </tr>
              <tr>
                <td>Dataset Description:</td>
                <td>{datasetData?.description ?? '(no description)'}</td>
              </tr>
              <tr>
                <td>Number of Reactions in Dataset:</td>
                <td>{datasetData?.num_reactions}</td>
              </tr>
            </tbody>
          </table>

          <div className="search-results-section">
            {error ? (
              <div className="no-results">
                <div className="title">Failed to load reactions: {error.message}</div>
              </div>
            ) : loading ? (
              <div className="loading">
                <LoadingSpinner />
              </div>
            ) : searchResults.length > 0 ? (
              <SearchResults
                searchResults={searchResults}
                isOverflow={(datasetData?.num_reactions ?? 0) > searchResults.length}
              />
            ) : (
              <div className="no-results">
                <div className="title">This dataset contains no reactions.</div>
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
};

export default MainDatasetView;
