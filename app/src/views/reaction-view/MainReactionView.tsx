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

// Read-only reaction page, structured like ord-app's ReactionPage: breadcrumbs,
// actions row, a header card with the reaction id + scheme preview, a
// Tabs/List toggle, and the section content.
import React, { useMemo, useState } from 'react';
import { useParams } from 'wouter';
import {
  Button,
  Flex,
  Paper,
  SegmentedControl,
  Title,
  type SegmentedControlItem,
} from '@mantine/core';
import { IconDownload } from '@tabler/icons-react';
import { useQuery } from '@tanstack/react-query';
import { ord } from 'ord-schema-protobufjs';
import Breadcrumbs, { type Breadcrumb } from '../../components/Breadcrumbs';
import CopyButton from '../../components/CopyButton';
import LoadingSpinner from '../../components/LoadingSpinner';
import { api } from '../../utils/api';
import { decodeReaction } from '../../utils/proto';
import { NotificationVariant, showNotification } from '../../utils/showNotification';
import { reactionViewContext } from './context';
import { ReactionContent } from './ReactionContent';
import ReactionPreviewStrip from './ReactionPreviewStrip';
import classes from './MainReactionView.module.scss';

const VIEW_MODE_OPTIONS: Array<SegmentedControlItem> = [
  { label: 'Tabs', value: 'tabs' },
  { label: 'List', value: 'list' },
];

interface ReactionRecord {
  reaction: ord.Reaction;
  datasetId?: string;
}

const MainReactionView: React.FC = () => {
  const { reactionId } = useParams<{ reactionId: string }>();
  const [viewMode, setViewMode] = useState<'tabs' | 'list'>('tabs');

  const { data, isLoading, error } = useQuery<ReactionRecord>({
    queryKey: ['reaction', reactionId],
    enabled: !!reactionId,
    queryFn: async () => {
      const response = await api.post<Array<{ proto: string; dataset_id?: string }>>(
        '/reactions',
        { reaction_ids: [reactionId] },
      );
      const record = response.data?.[0];
      if (!record?.proto) throw new Error(`Reaction ${reactionId} not found`);
      return {
        reaction: decodeReaction(record.proto),
        datasetId: record.dataset_id,
      };
    },
  });

  const breadcrumbs = useMemo((): Breadcrumb[] => {
    const crumbs: Breadcrumb[] = [{ title: 'Datasets', path: '/browse' }];
    if (data?.datasetId) {
      crumbs.push({ title: data.datasetId, path: `/dataset/${data.datasetId}` });
    }
    crumbs.push({ title: reactionId ?? '', path: `/id/${reactionId}` });
    return crumbs;
  }, [data?.datasetId, reactionId]);

  const downloadReaction = async () => {
    try {
      const response = await api.post<Blob>(
        '/download_search_results',
        { reaction_ids: [reactionId] },
        { responseType: 'blob' },
      );
      const url = URL.createObjectURL(response.data);
      const link = document.createElement('a');
      link.href = url;
      link.download = `${reactionId}.pb.gz`;
      link.click();
      setTimeout(() => {
        URL.revokeObjectURL(url);
        link.remove();
      }, 100);
    } catch (err) {
      console.error('Error downloading reaction:', err);
      showNotification({
        variant: NotificationVariant.ERROR,
        message: 'Download failed',
      });
    }
  };

  if (isLoading) {
    return <LoadingSpinner />;
  }

  if (error || !data) {
    return (
      <Paper
        radius="sm"
        p="xl"
        className={classes.emptyState}
      >
        Failed to load reaction {reactionId}.
      </Paper>
    );
  }

  return (
    <reactionViewContext.Provider value={{ reaction: data.reaction }}>
      <Breadcrumbs items={breadcrumbs} />
      <Flex
        direction="column"
        gap="sm"
      >
        <Flex
          justify="flex-end"
          align="flex-end"
        >
          <Button
            variant="transparent"
            leftSection={<IconDownload size={16} />}
            onClick={downloadReaction}
          >
            Download Reaction
          </Button>
        </Flex>

        <Paper
          radius="md"
          p="lg"
        >
          <Flex
            direction="column"
            gap="sm"
          >
            <Flex
              justify="space-between"
              align="center"
            >
              <Flex
                align="center"
                gap="sm"
                className={classes.titleWrapper}
              >
                <Title
                  className={classes.title}
                  order={2}
                >
                  {reactionId}
                </Title>
                <CopyButton textToCopy={reactionId ?? ''} />
              </Flex>
            </Flex>
            <ReactionPreviewStrip reaction={data.reaction} />
          </Flex>
        </Paper>

        <Paper className={classes.tableContainer}>
          <SegmentedControl
            value={viewMode}
            className={classes.controlBlock}
            color="var(--color-blue)"
            onChange={value => setViewMode(value as 'tabs' | 'list')}
            data={VIEW_MODE_OPTIONS}
          />
        </Paper>

        <Paper
          radius="md"
          p="lg"
        >
          <ReactionContent viewMode={viewMode} />
        </Paper>
      </Flex>
    </reactionViewContext.Provider>
  );
};

export default MainReactionView;
