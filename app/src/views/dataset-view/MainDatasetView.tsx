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

import React, { useState, useEffect } from 'react';
import { useParams } from 'react-router-dom';
import LoadingSpinner from '../../components/LoadingSpinner';
import ChartView from './ChartView';
import SearchResults from './SearchResults';
import { base64ToBytes } from '../../utils/base64';
import './MainDatasetView.scss';

interface Dataset {
  dataset_id: string;
  name?: string;
  description?: string;
  num_reactions: number;
}

interface SearchResult {
  reaction_id: string;
  proto: string;
  data?: any;
}

const MainDatasetView: React.FC = () => {
  const { datasetId } = useParams<{ datasetId: string }>();
  const [searchResults, setSearchResults] = useState<SearchResult[]>([]);
  const [loading, setLoading] = useState(true);
  const [datasetData, setDatasetData] = useState<Dataset | null>(null);
  const [isCollapsed, setIsCollapsed] = useState(true);
  const [searchTaskId, setSearchTaskId] = useState<string | null>(null);

  const expandOrShrink = () => {
    setIsCollapsed(!isCollapsed);
  };

  const getSearchResults = async () => {
    if (!datasetId) return;
    
    setLoading(true);
    try {
      // If this is the first time we are attempting the search, set the search task ID.
      let taskId = searchTaskId;
      if (taskId === null) {
        // Submit a search task to the server. This returns a GUID task ID.
        const taskRes = await fetch(`/api/submit_query?dataset_id=${datasetId}&limit=100`, {method: "GET"});
        taskId = await taskRes.json();
        setSearchTaskId(taskId);
      }
      
      // Check the status of the search task.
      const queryRes = await fetch(`/api/fetch_query_result?task_id=${taskId}`, {method: "GET"});
      
      if (queryRes.status === 200) {
        const searchResultsData = await queryRes.json();
        setSearchTaskId(null);
        
        // Deserialize the results and unpack protobuf for each reaction
        const processedResults = searchResultsData.map((reaction: SearchResult) => {
          // Convert base64 to bytes for protobuf parsing
          base64ToBytes(reaction.proto);
          // Note: This would require ord-schema package for full protobuf deserialization
          reaction.data = {
            identifiersList: [{ value: `Reaction ${reaction.reaction_id}` }]
          };
          return reaction;
        });
        
        setSearchResults(processedResults);
        setLoading(false);
      } else if (queryRes.status === 400) {
        setSearchTaskId(null);
        throw new Error(`Error - Search task ID ${taskId} does not exist`);
      } else if (queryRes.status === 500) {
        setSearchTaskId(null);
        throw new Error(`Error - Search task ID ${taskId} failed due to server error`);
      } else if (queryRes.status === 202) {
        // Search is still in progress, we'll poll again
        setTimeout(() => getSearchResults(), 1000);
      }
    } catch (error) {
      console.error('Search error:', error);
      setSearchResults([]);
      setLoading(false);
    }
  };

  // Fetch dataset metadata
  useEffect(() => {
    if (!datasetId) return;
    
    fetch(`/api/datasets`)
      .then(response => response.json())
      .then((datasets: Dataset[]) => {
        const dataset = datasets.find(d => d.dataset_id === datasetId);
        setDatasetData(dataset || null);
      })
      .catch(error => console.error('Error fetching dataset metadata:', error));
  }, [datasetId]);

  // Initial search results fetch
  useEffect(() => {
    if (datasetId) {
      getSearchResults();
    }
  }, [datasetId]);

  return (
    <div id="dataset-main">
      <h1>Dataset View</h1>
      <div 
        id="charts" 
        style={{
          display: isCollapsed ? 'grid' : 'flex',
          flexDirection: isCollapsed ? undefined : 'column'
        }}
      >
        <div 
          id="chartsection" 
          style={{ width: isCollapsed ? '40%' : '100%' }}
        >
          <div 
            id="chartsectioncharts" 
            style={{ display: isCollapsed ? 'block' : 'flex' }}
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
          <div id="expand">
            <button onClick={expandOrShrink}>
              <i 
                className="material-icons" 
                style={{ marginTop: '15%' }} 
                title={isCollapsed ? "Expand" : "Collapse"}
              >
                {isCollapsed ? "keyboard_double_arrow_right" : "keyboard_double_arrow_left"}
              </i>
            </button>
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
            {!loading && searchResults?.length > 0 ? (
              <SearchResults 
                searchResults={searchResults}
                isOverflow={false}
              />
            ) : !loading && !searchResults?.length ? (
              <div className="no-results">
                <div className="title">This dataset contains no reactions.</div>
              </div>
            ) : (
              <div className="loading">
                <LoadingSpinner />
              </div>
            )}
          </div>
        </div>
      </div>
    </div>
  );
};

export default MainDatasetView;