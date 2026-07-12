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

import { useEffect, useMemo, useState, type ReactNode } from 'react';
import { Badge, TextInput, Title } from '@mantine/core';
import { IconSearch } from '@tabler/icons-react';
import Pagination from './Pagination';
import classes from './EntityTable.module.scss';

interface EntityTableProps<T> {
  tableData: T[];
  title?: string;
  displaySearch?: boolean;
  children: (entities: T[]) => ReactNode;
}

/**
 * Generic client-side filtered + paginated list wrapper. The render prop
 * receives the current page of entities; header shows a title, count badge,
 * and free-text filter in the ord-app style.
 */
function EntityTable<T>({
  tableData,
  title = '',
  displaySearch = true,
  children,
}: EntityTableProps<T>) {
  const [searchString, setSearchString] = useState('');
  const [currentPage, setCurrentPage] = useState(1);
  const [pageSize, setPageSize] = useState(10);

  // Back to the first page whenever the page no longer exists.
  useEffect(() => {
    setCurrentPage(1);
  }, [pageSize, searchString, tableData]);

  const searchArray = useMemo(() => {
    return searchString.split(' ').filter(term => term.trim() !== '');
  }, [searchString]);

  const filteredEntities = useMemo(() => {
    if (!searchArray.length) return tableData;
    return tableData.filter(entity => {
      const fields = entity as Record<string, unknown>;
      return searchArray.every(param =>
        Object.values(fields).some(value => {
          const stringified = value == null ? 'null' : String(value);
          return stringified.toLowerCase().includes(param.toLowerCase());
        }),
      );
    });
  }, [tableData, searchArray]);

  const lastPage = useMemo(
    () => Math.max(1, Math.ceil(filteredEntities.length / pageSize)),
    [filteredEntities, pageSize],
  );

  const paginatedEntities = useMemo(() => {
    const start = (currentPage - 1) * pageSize;
    return filteredEntities.slice(start, start + pageSize);
  }, [filteredEntities, currentPage, pageSize]);

  if (!tableData.length) {
    return null;
  }

  return (
    <div className={classes.tableMain}>
      {(title || displaySearch) && (
        <div className={classes.header}>
          {title && (
            <div className={classes.titleHolder}>
              <Title order={2}>{title}</Title>
              <Badge
                className={classes.countBadge}
                variant="light"
                color="gray"
                radius="xl"
              >
                {filteredEntities.length.toLocaleString()}
              </Badge>
            </div>
          )}
          {displaySearch && (
            <TextInput
              classNames={{ root: classes.search, input: classes.searchInput }}
              leftSection={<IconSearch size={14} />}
              placeholder="Filter…"
              value={searchString}
              onChange={event => setSearchString(event.currentTarget.value)}
              aria-label="Filter table"
            />
          )}
        </div>
      )}

      <div className={classes.content}>{children(paginatedEntities)}</div>

      {lastPage > 1 || filteredEntities.length > 10 ? (
        <Pagination
          currentPage={currentPage}
          totalPages={lastPage}
          onPageChange={setCurrentPage}
          rowsPerPage={pageSize}
          onRowsPerPageChange={setPageSize}
        />
      ) : null}
    </div>
  );
}

export default EntityTable;
