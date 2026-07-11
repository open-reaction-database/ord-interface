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

import React from 'react';
import { Button, Group, Pagination as MantinePagination, Select } from '@mantine/core';
import { IconChevronLeft, IconChevronRight } from '@tabler/icons-react';
import classes from './Pagination.module.scss';

interface PaginationProps {
  currentPage: number;
  totalPages: number;
  onPageChange: (page: number) => void;
  rowsPerPage: number;
  onRowsPerPageChange: (rows: number) => void;
}

const ROW_PER_PAGE_OPTIONS = ['10', '25', '50', '100'];

/** Pager matching ord-app's house Pagination component. */
const Pagination: React.FC<PaginationProps> = ({
  currentPage,
  totalPages,
  onPageChange,
  rowsPerPage,
  onRowsPerPageChange,
}) => {
  return (
    <Group
      justify="space-between"
      className={classes.paginationContainer}
    >
      <Group align="center">
        <Button
          className={classes.controlButton}
          leftSection={<IconChevronLeft size={16} />}
          variant="transparent"
          disabled={currentPage === 1}
          onClick={() => onPageChange(currentPage - 1)}
        >
          Previous
        </Button>

        <MantinePagination
          className={classes.paginationRoot}
          total={totalPages}
          value={currentPage}
          onChange={onPageChange}
          size="lg"
          withControls={false}
          siblings={0}
        />

        <Button
          className={classes.controlButton}
          rightSection={<IconChevronRight size={16} />}
          variant="transparent"
          disabled={currentPage === totalPages}
          onClick={() => onPageChange(currentPage + 1)}
        >
          Next
        </Button>
      </Group>

      <Select
        classNames={{ root: classes.select, input: classes.selectInput }}
        value={String(rowsPerPage)}
        onChange={value => onRowsPerPageChange(Number.parseInt(value as string))}
        allowDeselect={false}
        data={ROW_PER_PAGE_OPTIONS}
        rightSectionWidth={30}
      />
    </Group>
  );
};

export default Pagination;
