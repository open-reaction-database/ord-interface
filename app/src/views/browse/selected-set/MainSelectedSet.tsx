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
import { useSearch } from 'wouter';
import { Badge, Button, Paper, Text, Title } from '@mantine/core';
import { IconDownload } from '@tabler/icons-react';
import LoadingSpinner from '../../../components/LoadingSpinner';
import ReactionCard from '../../../components/ReactionCard';
import DownloadResults from '../../../components/DownloadResults';
import CopyButton from '../../../components/CopyButton';
import { api } from '../../../utils/api';
import { decodeReaction } from '../../../utils/proto';
import type { SearchResult } from '../../../types/search';
import classes from './MainSelectedSet.module.scss';

const MainSelectedSet: React.FC = () => {
  const search = useSearch();
  const [reactions, setReactions] = useState<SearchResult[]>([]);
  const [loading, setLoading] = useState(true);
  const [showDownloadResults, setShowDownloadResults] = useState(false);

  const reactionIds = useMemo(
    () => new URLSearchParams(search).getAll('reaction_id'),
    [search],
  );
  const reactionIdsKey = reactionIds.join(',');
  const fullUrl = window.location.href;

  const getSelectedReactions = useCallback(async () => {
    setLoading(true);
    try {
      const response = await api.post<Array<Omit<SearchResult, 'data'>>>('/reactions', {
        reaction_ids: reactionIds,
      });

      const decoded: SearchResult[] = response.data.map(r => ({
        ...r,
        data: decodeReaction(r.proto),
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
    return <LoadingSpinner />;
  }

  if (!reactions.length) {
    return (
      <Paper
        radius="sm"
        p="xl"
        className={classes.emptyState}
      >
        <Text c="secondary.1">
          There was an issue fetching your selected reactions.
        </Text>
      </Paper>
    );
  }

  return (
    <div className={classes.selectedSetMain}>
      <div className={classes.header}>
        <div className={classes.titleHolder}>
          <Title order={2}>Reaction set</Title>
          <Badge
            className={classes.countBadge}
            variant="light"
            color="gray"
            radius="xl"
          >
            {reactions.length.toLocaleString()}
          </Badge>
        </div>
        <div className={classes.actions}>
          <CopyButton
            textToCopy={fullUrl}
            icon="share"
            buttonText="Shareable link"
          />
          <Button
            variant="default"
            radius="sm"
            leftSection={<IconDownload size={16} />}
            disabled={!reactionIds.length}
            onClick={() => setShowDownloadResults(true)}
          >
            Download reaction set
          </Button>
        </div>
      </div>

      <div>
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
