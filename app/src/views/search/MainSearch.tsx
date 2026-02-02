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

import React, { useState, useEffect, useCallback, useRef, useMemo } from 'react';
import { useLocation, useNavigate } from 'react-router-dom';
import SearchOptions from './SearchOptions';
import SearchResults from './SearchResults';
import LoadingSpinner from '../../components/LoadingSpinner';
import reaction_pb from 'ord-schema';
import { base64ToBytes } from '../../utils/base64';
import './MainSearch.scss';

interface SearchResult {
  reaction_id: string;
  proto: string;
  data: any;
  [key: string]: any;
}

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

const MainSearch: React.FC = () => {
  const location = useLocation();
  const navigate = useNavigate();
  const [searchResults, setSearchResults] = useState<SearchResult[]>([]);
  const [loading, setLoading] = useState(false);
  const [showOptions, setShowOptions] = useState(false);
  const [searchLoadStatus, setSearchLoadStatus] = useState<Response | null>(null);
  const [searchPollingInterval, setSearchPollingInterval] = useState<number | null>(null);
  const [searchTaskId, setSearchTaskId] = useState<string | null>(null);
  const searchPollingIntervalRef = useRef<number | null>(null);
  const searchTaskIdRef = useRef<string | null>(null);

  // Check if there are any search parameters
  const hasSearchParams = useMemo(() => {
    const params = new URLSearchParams(location.search);
    // Check for any meaningful search parameters (excluding just 'limit')
    const meaningfulParams = ['component', 'dataset_id', 'doi', 'reaction_id', 'reaction_smarts', 
                              'min_yield', 'max_yield', 'min_conversion', 'max_conversion'];
    return meaningfulParams.some(param => params.has(param));
  }, [location.search]);

  const getSearchResults = useCallback(async () => {
    setLoading(true);
    const urlQuery = location.search;
    
    try {
      let currentTaskId = searchTaskIdRef.current;
      
      // If this is the first time we are attempting the search, set the search task ID.
      if (currentTaskId == null) {
        // Submit a search task to the server. This returns a GUID task ID.
        const taskres = await fetch(`/api/submit_query${urlQuery}`, { method: 'GET' });
        const taskId = await taskres.json();
        setSearchTaskId(taskId);
        searchTaskIdRef.current = taskId;
        currentTaskId = taskId;
      }
      
      // Check the status of the search task.
      const queryRes = await fetch(`/api/fetch_query_result?task_id=${currentTaskId}`, { method: 'GET' });
      setSearchLoadStatus(queryRes);
      
      // If one of these codes, search is finished. Return the results or lack of results.
      if (queryRes?.status === 200) {
        setSearchTaskId(null);
        searchTaskIdRef.current = null;
        if (searchPollingIntervalRef.current) {
          clearInterval(searchPollingIntervalRef.current);
          setSearchPollingInterval(null);
          searchPollingIntervalRef.current = null;
        }
        
        const searchResultsData = await queryRes.json();
        
        // Deserialize the results.
        const results: SearchResult[] = searchResultsData;
        // unpack protobuff for each reaction in results
        results.forEach((reaction) => {
          const bytes = base64ToBytes(reaction.proto);
          reaction.data = reaction_pb.Reaction.deserializeBinary(new Uint8Array(bytes)).toObject();
        });
        
        setSearchResults(results);
        setLoading(false);
        return queryRes.status;
      } else if (queryRes?.status === 404) {
        setSearchTaskId(null);
        searchTaskIdRef.current = null;
        if (searchPollingIntervalRef.current) {
          clearInterval(searchPollingIntervalRef.current);
          setSearchPollingInterval(null);
          searchPollingIntervalRef.current = null;
        }
        throw new Error(`Error - Search task ID ${currentTaskId} does not exist`);
      } else if (queryRes?.status >= 500) {
        setSearchTaskId(null);
        searchTaskIdRef.current = null;
        if (searchPollingIntervalRef.current) {
          clearInterval(searchPollingIntervalRef.current);
          setSearchPollingInterval(null);
          searchPollingIntervalRef.current = null;
        }
        throw new Error(`Error - Search task ID ${currentTaskId} failed due to server error`);
      }
      
      return queryRes?.status;
    } catch (e) {
      console.log(e);
      setSearchResults([]);
      setLoading(false);
      return null;
    }
  }, [location.search]);

  const updateSearchOptions = (options: SearchOptionsData) => {
    const searchParams = new URLSearchParams();

    // reagent options
    if (options.reagent.reagents.length) {
      options.reagent.reagents.forEach(reagent => {
        searchParams.append('component', `${reagent.smileSmart};${reagent.source};${reagent.matchMode.toLowerCase()}`);
      });

      searchParams.set('use_stereochemistry', options.reagent.useStereochemistry.toString());
      searchParams.set('similarity', options.reagent.similarityThreshold.toString());
    }

    // dataset options
    if (options.dataset.datasetIds.length) {
      options.dataset.datasetIds.forEach(id => {
        searchParams.append('dataset_id', id);
      });
    }
    if (options.dataset.DOIs.length) {
      options.dataset.DOIs.forEach(doi => {
        searchParams.append('doi', doi);
      });
    }

    // reaction options
    if (options.reaction.reactionIds.length) {
      options.reaction.reactionIds.forEach(id => {
        searchParams.append('reaction_id', id);
      });
    }
    if (options.reaction.reactionSmarts.length) {
      options.reaction.reactionSmarts.forEach(smarts => {
        searchParams.append('reaction_smarts', smarts);
      });
    }
    
    // yield and conversion add if not max values, otherwise remove from query
    if (options.reaction.min_yield !== 0 || options.reaction.max_yield !== 100) {
      searchParams.set('min_yield', options.reaction.min_yield.toString());
      searchParams.set('max_yield', options.reaction.max_yield.toString());
    }
    if (options.reaction.min_conversion !== 0 || options.reaction.max_conversion !== 100) {
      searchParams.set('min_conversion', options.reaction.min_conversion.toString());
      searchParams.set('max_conversion', options.reaction.max_conversion.toString());
    }

    // general options
    searchParams.set('limit', options.general.limit.toString() || '100');

    // navigate to search page with new params
    navigate({ pathname: '/search', search: searchParams.toString() });
  };

  useEffect(() => {
    // Clear any existing task and interval when search parameters change
    setSearchTaskId(null);
    searchTaskIdRef.current = null;
    setSearchResults([]);
    if (searchPollingIntervalRef.current) {
      clearInterval(searchPollingIntervalRef.current);
      setSearchPollingInterval(null);
      searchPollingIntervalRef.current = null;
    }

    // Only fetch results if there are search parameters
    if (!hasSearchParams) {
      setLoading(false);
      return;
    }

    const fetchResults = async () => {
      // Fetch results. If server returns a 202, set up a poll to keep checking back until we have results.
      const status = await getSearchResults();
      
      // If status is 202, set up polling
      if (status === 202) {
        // Clear any existing interval first
        if (searchPollingIntervalRef.current) {
          clearInterval(searchPollingIntervalRef.current);
        }
        
        const interval = setInterval(() => {
          getSearchResults();
        }, 1000);
        
        setSearchPollingInterval(interval);
        searchPollingIntervalRef.current = interval;
        
        // Set timeout to stop polling after 2 minutes
        setTimeout(() => {
          if (searchPollingIntervalRef.current === interval) {
            clearInterval(interval);
            setSearchTaskId(null);
            searchTaskIdRef.current = null;
            setLoading(false);
            setSearchPollingInterval(null);
            searchPollingIntervalRef.current = null;
          }
        }, 120000);
      }
    };

    fetchResults();
  }, [location.search, getSearchResults, hasSearchParams]); // Re-run when search parameters change

  // Cleanup interval on unmount
  useEffect(() => {
    return () => {
      if (searchPollingIntervalRef.current) {
        clearInterval(searchPollingIntervalRef.current);
      }
    };
  }, []);

  return (
    <div id="search-main">
      <div className={`search-options-container ${showOptions ? 'slide-out' : 'hidden'}`}>
        <div className="title">Filters & Options</div>
        <div className="options-holder">
          <SearchOptions onSearchOptions={updateSearchOptions} />
        </div>
        <div className="slide-out-tab" onClick={() => setShowOptions(!showOptions)}>
          <div className="line"></div>
          <div className="line"></div>
          <div className="line"></div>
        </div>
      </div>
      
      <div className="search-results">
        {!hasSearchParams && (
          <div className="no-results">
            <div className="title">Enter search criteria using the filters and options panel, then click Search.</div>
          </div>
        )}
        {hasSearchParams && !loading && searchResults?.length > 0 && (
          <SearchResults searchResults={searchResults} />
        )}
        {hasSearchParams && !loading && !searchResults?.length && (
          <div className="no-results">
            <div className="title">No results. Adjust the filters and options and search again.</div>
          </div>
        )}
        {hasSearchParams && loading && (
          <div className="loading">
            <LoadingSpinner />
          </div>
        )}
      </div>
    </div>
  );
};

export default MainSearch;