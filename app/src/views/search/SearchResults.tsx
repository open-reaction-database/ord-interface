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

import React, { useEffect, useState } from 'react';
import { useLocation, useSearch } from 'wouter';
import { Button } from '@mantine/core';
import { IconDownload, IconStack2 } from '@tabler/icons-react';
import EntityTable from '../../components/EntityTable';
import ReactionCard from '../../components/ReactionCard';
import CopyButton from '../../components/CopyButton';
import DownloadResults from '../../components/DownloadResults';
import type { SearchResult } from '../../types/search';
import classes from './SearchResults.module.scss';

interface SearchResultsProps {
  searchResults: SearchResult[];
}

const SearchResults: React.FC<SearchResultsProps> = ({ searchResults }) => {
  const search = useSearch();
  const [, navigate] = useLocation();
  const [selectedReactions, setSelectedReactions] = useState<string[]>([]);
  const [showDownloadResults, setShowDownloadResults] = useState(false);

  const locationSearch = search ? `?${search}` : '';

  const updateSelectedReactions = (reactionId: string, isSelected: boolean) => {
    if (isSelected) {
      setSelectedReactions(prev => [...prev, reactionId]);
    } else {
      setSelectedReactions(prev => prev.filter(id => id !== reactionId));
    }
  };

  const goToViewSelected = () => {
    // Persist selection so it survives navigating back from the selected-set page.
    localStorage.setItem(
      'storedSet',
      JSON.stringify({
        query: locationSearch,
        reactions: selectedReactions,
      }),
    );

    const params = new URLSearchParams();
    selectedReactions.forEach(id => params.append('reaction_id', id));
    navigate(`/selected-set?${params.toString()}`);
  };

  // Restore the prior selection when the user returns from /selected-set to the
  // same search URL. Driven by the query string only; results render directly
  // from the prop to avoid a one-frame empty flash on mount.
  useEffect(() => {
    const storedSetStr = localStorage.getItem('storedSet');
    if (!storedSetStr) return;
    try {
      const storedSet = JSON.parse(storedSetStr);
      if (locationSearch === storedSet.query) {
        setSelectedReactions(storedSet.reactions || []);
      }
    } catch (error) {
      console.error('Error parsing stored set:', error);
    }
  }, [locationSearch]);

  return (
    <div className={classes.searchResultsMain}>
      {searchResults.length > 0 && (
        <EntityTable
          tableData={searchResults}
          title="Search results"
          displaySearch={false}
        >
          {entities => (
            <>
              <div className={classes.actionButtons}>
                <CopyButton
                  textToCopy={window.location.href}
                  icon="share"
                  buttonText="Shareable link"
                />
                <Button
                  variant="default"
                  radius="sm"
                  leftSection={<IconDownload size={16} />}
                  disabled={!searchResults.length}
                  onClick={() => setShowDownloadResults(true)}
                >
                  Download all results
                </Button>
              </div>
              {entities.map(row => (
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
        <div className={classes.viewSelectedContainer}>
          <Button
            radius="xl"
            leftSection={<IconStack2 size={16} />}
            onClick={goToViewSelected}
          >
            View {selectedReactions.length} selected reaction
            {selectedReactions.length === 1 ? '' : 's'}
          </Button>
        </div>
      )}

      <DownloadResults
        reactionIds={searchResults.map(result => result.reaction_id)}
        showDownloadResults={showDownloadResults}
        onHideDownloadResults={() => setShowDownloadResults(false)}
      />
    </div>
  );
};

export default SearchResults;
