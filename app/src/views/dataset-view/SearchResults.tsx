import React, { useState, useEffect } from 'react';
import { useNavigate } from 'react-router-dom';
import EntityTable from '../../components/EntityTable';
import ReactionCard from '../../components/ReactionCard';
import CopyButton from '../../components/CopyButton';
import DownloadResults from '../../components/DownloadResults';
import './SearchResults.scss';

interface SearchResult {
  reaction_id: string;
  [key: string]: any;
}

interface SearchResultsProps {
  searchResults: SearchResult[];
  isOverflow: boolean;
}

const SearchResults: React.FC<SearchResultsProps> = ({
  searchResults,
  isOverflow
}) => {
  const navigate = useNavigate();
  const [formattedResults, setFormattedResults] = useState<SearchResult[]>([]);
  const [selectedReactions, setSelectedReactions] = useState<string[]>([]);
  const [showDownloadResults, setShowDownloadResults] = useState(false);



  const goToViewSelected = () => {
    // Store selected reactions in sessionStorage for navigation
    sessionStorage.setItem('selectedReactions', JSON.stringify({
      query: window.location.search,
      reactions: selectedReactions
    }));
    
    navigate(`/selected-set?reaction_ids=${selectedReactions.join(',')}`);
  };

  useEffect(() => {
    setFormattedResults(searchResults);
    
    // Check if there are stored selected reactions for this query
    const stored = sessionStorage.getItem('selectedReactions');
    if (stored) {
      const { query, reactions } = JSON.parse(stored);
      if (query === window.location.search) {
        setSelectedReactions(reactions);
      }
    }
  }, [searchResults]);

  const title = isOverflow 
    ? "100 Reactions From This Dataset (Sample)" 
    : `Reactions in this Dataset (${formattedResults.length} Reactions)`;

  return (
    <div className="search-results">
      {formattedResults.length > 0 && (
        <EntityTable
          tableData={formattedResults}
          title={title}
          displaySearch={false}
        >
          {(entities) => (
            <>
              <div className="search-results__action-buttons">
                <CopyButton
                  textToCopy={window.location.href}
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
              
              {entities.map((row: SearchResult) => (
                <ReactionCard
                  key={row.reaction_id}
                  reaction={row}
                  isSelectable={false}
                  isSelected={selectedReactions.includes(row.reaction_id)}
                  onSelectionChange={(reactionId, isSelected) => {
                    if (isSelected) {
                      setSelectedReactions(prev => [...prev, reactionId]);
                    } else {
                      setSelectedReactions(prev => prev.filter(id => id !== reactionId));
                    }
                  }}
                />
              ))}
            </>
          )}
        </EntityTable>
      )}

      {selectedReactions.length > 0 && (
        <div className="search-results__view-selected">
          <div 
            className="search-results__view-selected-button"
            onClick={goToViewSelected}
          >
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