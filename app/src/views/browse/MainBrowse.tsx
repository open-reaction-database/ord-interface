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

import React, { useMemo, useState } from 'react';
import { Link } from 'wouter';
import { ActionIcon, Badge, Paper, Text, Tooltip } from '@mantine/core';
import { IconArrowsDiagonal } from '@tabler/icons-react';
import type { MRT_ColumnDef } from 'mantine-react-table';
import { useQuery } from '@tanstack/react-query';
import { DataTable } from '../../components/DataTable';
import EntityTable from '../../components/EntityTable';
import FloatingModal from '../../components/FloatingModal';
import LoadingSpinner from '../../components/LoadingSpinner';
import { api } from '../../utils/api';
import type { Dataset } from '../../types/search';
import classes from './MainBrowse.module.scss';

// Clamped description with an expand button that opens the dataset's full
// details in a modal, so the table layout never reflows.
const DescriptionCell: React.FC<{ dataset: Dataset }> = ({ dataset }) => {
  const [expanded, setExpanded] = useState(false);

  return (
    <div className={classes.descriptionCell}>
      <Text
        lineClamp={2}
        className={classes.descriptionText}
      >
        {dataset.description}
      </Text>
      <Tooltip label="Show dataset details">
        <ActionIcon
          variant="transparent"
          color="secondary.1"
          onClick={() => setExpanded(true)}
          aria-label="Show dataset details"
        >
          <IconArrowsDiagonal size={16} />
        </ActionIcon>
      </Tooltip>
      {expanded && (
        <FloatingModal
          title={dataset.name || dataset.dataset_id}
          size="lg"
          onCloseModal={() => setExpanded(false)}
        >
          <div className={classes.modalDatasetId}>{dataset.dataset_id}</div>
          <div className={classes.modalMeta}>
            <span>
              <strong>Submitted:</strong> {dataset.submitted_at ?? '—'}
            </span>
            <span>
              <strong>Reactions:</strong> {dataset.num_reactions.toLocaleString()}
            </span>
          </div>
          <div className={classes.modalDescription}>{dataset.description}</div>
        </FloatingModal>
      )}
    </div>
  );
};

const MainBrowse: React.FC = () => {
  const {
    data: tableData,
    isLoading,
    error,
  } = useQuery<Dataset[]>({
    queryKey: ['datasets'],
    queryFn: () => api.get<Dataset[]>('/datasets').then(res => res.data),
  });

  const columns = useMemo<Array<MRT_ColumnDef<Dataset>>>(
    () => [
      {
        accessorKey: 'dataset_id',
        header: 'Dataset ID',
        size: 300,
        Cell: ({ row }) => (
          <Link
            href={`/dataset/${row.original.dataset_id}`}
            className={classes.datasetLink}
          >
            {row.original.dataset_id}
          </Link>
        ),
      },
      {
        accessorKey: 'name',
        header: 'Name',
        size: 260,
      },
      {
        accessorKey: 'description',
        header: 'Description',
        enableSorting: false,
        Cell: ({ row }) => <DescriptionCell dataset={row.original} />,
      },
      {
        accessorKey: 'num_reactions',
        header: 'Size',
        size: 120,
        Cell: ({ row }) => (
          <Badge
            className={classes.sizeBadge}
            variant="light"
            color="gray"
            radius="xl"
          >
            {row.original.num_reactions.toLocaleString()}
          </Badge>
        ),
      },
    ],
    [],
  );

  if (isLoading) {
    return <LoadingSpinner />;
  }

  if (error || !tableData?.length) {
    return (
      <Paper
        radius="sm"
        p="xl"
        className={classes.errorState}
      >
        <Text c="secondary.1">
          {error ? `Failed to load datasets: ${error.message}` : 'No datasets found.'}
        </Text>
      </Paper>
    );
  }

  return (
    <Paper
      radius="sm"
      p="lg"
      className={classes.browsePaper}
    >
      <EntityTable
        tableData={tableData}
        title="Datasets"
      >
        {(entities: Dataset[]) => (
          <DataTable
            columns={columns}
            data={entities}
            enableSorting={false}
          />
        )}
      </EntityTable>
    </Paper>
  );
};

export default MainBrowse;
