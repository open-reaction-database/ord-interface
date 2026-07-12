/**
 * Copyright 2026 Open Reaction Database Project Authors
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
import { Button } from '@mantine/core';
import { IconDownload } from '@tabler/icons-react';
import EntityTable from '../../components/EntityTable';
import ReactionCard from '../../components/ReactionCard';
import CopyButton from '../../components/CopyButton';
import DownloadResults from '../../components/DownloadResults';
import type { SearchResult } from '../../types/search';
import classes from './SearchResults.module.scss';

interface SearchResultsProps {
  searchResults: SearchResult[];
  isOverflow: boolean;
}

const SearchResults: React.FC<SearchResultsProps> = ({ searchResults, isOverflow }) => {
  const [showDownloadResults, setShowDownloadResults] = useState(false);

  const title = isOverflow ? 'Reactions (first 100 shown)' : 'Reactions';

  if (searchResults.length === 0) {
    return null;
  }

  return (
    <div className={classes.searchResults}>
      <EntityTable
        tableData={searchResults}
        title={title}
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
                isSelectable={false}
              />
            ))}
          </>
        )}
      </EntityTable>

      <DownloadResults
        reactionIds={searchResults.map(result => result.reaction_id)}
        showDownloadResults={showDownloadResults}
        onHideDownloadResults={() => setShowDownloadResults(false)}
      />
    </div>
  );
};

export default SearchResults;
