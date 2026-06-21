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

import React, { useState, useEffect, useMemo, useCallback } from 'react';
import { useSearchParams } from 'react-router-dom';
import reaction_pb from 'ord-schema';
import LoadingSpinner from '../../../components/LoadingSpinner';
import ReactionCard from '../../../components/ReactionCard';
import DownloadResults from '../../../components/DownloadResults';
import CopyButton from '../../../components/CopyButton';
import { base64ToBytes } from '../../../utils/base64';
import type { SearchResult } from '../../../types/search';
import './MainSelectedSet.scss';

const MainSelectedSet: React.FC = () => {
  const [searchParams] = useSearchParams();
  const [reactions, setReactions] = useState<SearchResult[]>([]);
  const [loading, setLoading] = useState(true);
  const [showDownloadResults, setShowDownloadResults] = useState(false);

  const reactionIds = useMemo(() => searchParams.getAll('reaction_id'), [searchParams]);
  const reactionIdsKey = reactionIds.join(',');
  const fullUrl = window.location.href;

  const getSelectedReactions = useCallback(async () => {
    setLoading(true);
    try {
      const response = await fetch('/api/reactions', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ reaction_ids: reactionIds }),
      });

      if (!response.ok) throw new Error('Failed to fetch reactions');
      const fetched = (await response.json()) as Array<Omit<SearchResult, 'data'>>;

      const decoded: SearchResult[] = fetched.map(r => ({
        ...r,
        data: reaction_pb.Reaction.deserializeBinary(
          new Uint8Array(base64ToBytes(r.proto)),
        ).toObject(),
      }));

      setReactions(decoded);
    } catch (error) {
      console.error('Error fetching reactions:', error);
      setReactions([]);
    } finally {
      setLoading(false);
    }
  }, [reactionIdsKey]); // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    if (reactionIds.length > 0) {
      getSelectedReactions();
    } else {
      setLoading(false);
    }
  }, [reactionIdsKey, getSelectedReactions]); // eslint-disable-line react-hooks/exhaustive-deps

  if (loading) {
    return (
      <div id="selected-set-main">
        <div className="loading">
          <LoadingSpinner />
        </div>
      </div>
    );
  }

  if (!reactions.length) {
    return (
      <div id="selected-set-main">
        <div className="header">
          <div className="title">Reaction Set</div>
        </div>
        <div className="no-results">
          <div className="title">
            There was an issue fetching your selected reactions.
          </div>
        </div>
      </div>
    );
  }

  return (
    <div id="selected-set-main">
      <div className="header">
        <div className="title">Reaction Set</div>
        <div className="action-button-holder">
          <CopyButton
            textToCopy={fullUrl}
            icon="share"
            buttonText="Shareable Link"
          />
          <button
            disabled={!reactionIds.length}
            onClick={() => setShowDownloadResults(true)}
          >
            Download Reaction Set
          </button>
        </div>
      </div>

      <div className="selected-set">
        {reactions.map(reaction => (
          <ReactionCard
            key={reaction.reaction_id}
            reaction={reaction}
            isSelectable={false}
            isSelected={false}
          />
        ))}
      </div>

      <DownloadResults
        reactionIds={reactionIds}
        showDownloadResults={showDownloadResults}
        onHideDownloadResults={() => setShowDownloadResults(false)}
      />
    </div>
  );
};

export default MainSelectedSet;
