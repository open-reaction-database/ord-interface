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

import {
  MantineReactTable,
  useMantineReactTable,
  type MRT_RowData,
  type MRT_TableOptions,
} from 'mantine-react-table';
import { IconChevronDown, IconChevronUp, IconSelector } from '@tabler/icons-react';
import classes from './DataTable.module.scss';

/** Table wrapper ported from ord-app's DataTable so tables look identical. */
export function DataTable<T extends MRT_RowData>({
  columns,
  data,
  ...rest
}: MRT_TableOptions<T>) {
  const table = useMantineReactTable<T>({
    columns,
    layoutMode: 'semantic',
    data: data,
    enableColumnActions: false,
    enableTopToolbar: false,
    enableBottomToolbar: false,
    enablePagination: false,
    enableStickyHeader: true,
    mantinePaperProps: {
      className: classes.paper,
      withBorder: false,
      ...rest.mantinePaperProps,
    },
    mantineTableContainerProps: {
      className: classes.tableContainer,
      ...rest.mantineTableContainerProps,
    },
    mantineTableProps: {
      className: classes.table,
      ...rest.mantineTableProps,
    },
    mantineTableHeadCellProps: {
      className: classes.headerCell,
      ...rest.mantineTableHeadCellProps,
    },
    mantineTableBodyCellProps: {
      className: classes.bodyCell,
      ...rest.mantineTableBodyCellProps,
    },
    mantineSkeletonProps: {
      height: 21.7,
      width: '90%',
      ...rest.mantineSkeletonProps,
    },
    icons: {
      IconSortDescending: IconChevronDown,
      IconSortAscending: IconChevronUp,
      IconArrowsSort: IconSelector,
      ...rest.icons,
    },
    ...rest,
  });

  return <MantineReactTable table={table} />;
}
