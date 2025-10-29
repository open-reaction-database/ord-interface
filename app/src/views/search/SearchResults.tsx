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

import React, { useState, useEffect, useMemo } from 'react';
import { useNavigate, useLocation } from 'react-router-dom';
import EntityTable from '../../components/EntityTable';
import ReactionCard from '../../components/ReactionCard';
import CopyButton from '../../components/CopyButton';
import DownloadResults from '../../components/DownloadResults';
import './SearchResults.scss';

interface SearchResult {
  reaction_id: string;
  proto: string;
  data: any;
  [key: string]: any;
}

interface SearchResultsProps {
  searchResults: SearchResult[];
}

const SearchResults: React.FC<SearchResultsProps> = ({ searchResults }) => {
  const navigate = useNavigate();
  const location = useLocation();
  const [formattedResults, setFormattedResults] = useState<SearchResult[]>([]);
  const [selectedReactions, setSelectedReactions] = useState<string[]>([]);
  const [showDownloadResults, setShowDownloadResults] = useState(false);

  const fullUrl = useMemo(() => {
    return window.location.href;
  }, []);

  const updateSelectedReactions = (reactionId: string, isSelected: boolean) => {
    if (isSelected) {
      setSelectedReactions(prev => [...prev, reactionId]);
    } else {
      setSelectedReactions(prev => prev.filter(id => id !== reactionId));
    }
  };

  const goToViewSelected = () => {
    // TODO: Implement Vuex store equivalent or local storage for storedSet
    // Store storedSet so we can retrieve it if user comes back from selected-set
    localStorage.setItem('storedSet', JSON.stringify({
      query: location.search,
      reactions: selectedReactions
    }));
    
    const params = new URLSearchParams();
    selectedReactions.forEach(id => params.append('reaction_id', id));
    navigate(`/selected-set?${params.toString()}`);
  };

  useEffect(() => {
    setFormattedResults(searchResults);
    
    // If query matches storedSet, set selectedReactions
    const storedSetStr = localStorage.getItem('storedSet');
    if (storedSetStr) {
      try {
        const storedSet = JSON.parse(storedSetStr);
        if (location.search === storedSet.query) {
          setSelectedReactions(storedSet.reactions || []);
        }
      } catch (error) {
        console.error('Error parsing stored set:', error);
      }
    }
  }, [searchResults, location.search]);

  return (
    <div className="search-results-main">
      {formattedResults.length > 0 && (
        <EntityTable
          tableData={formattedResults}
          title="Search Results"
          displaySearch={false}
        >
          {(entities: any[]) => (
            <>
              <div className="action-button-holder">
                <CopyButton
                  textToCopy={fullUrl}
                  icon="share"
                  buttonText="Shareable Link"
                />
                <button
                  disabled={!formattedResults.length}
                  onClick={() => setShowDownloadResults(true)}
                >
                  Download All Search Results
                </button>
              </div>
              {entities.map((row: any) => (
                <ReactionCard
                  key={row.reaction_id}
                  reaction={row}
                  isSelectable={true}
                  isSelected={selectedReactions.includes(row.reaction_id)}
                  onSelectionChange={updateSelectedReactions}
                />
              ))}
            </>
          )}
        </EntityTable>
      )}
      
      {selectedReactions.length > 0 && (
        <div className="view-selected-container">
          <div className="view-selected-button" onClick={goToViewSelected}>
            View {selectedReactions.length} selected reactions
          </div>
        </div>
      )}
      
      <DownloadResults
        reactionIds={formattedResults.map(result => result.reaction_id)}
        showDownloadResults={showDownloadResults}
        onHideDownloadResults={() => setShowDownloadResults(false)}
      />
    </div>
  );
};

export default SearchResults;