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

import React, { type ReactNode } from 'react';
import './EntityTable.scss';

interface EntityTableProps {
  tableData: any[];
  title?: string;
  children: (entities: any[]) => ReactNode;
}

const EntityTable: React.FC<EntityTableProps> = ({ 
  tableData, 
  title = '', 
  children 
}) => {
  if (!tableData.length) {
    return null;
  }

  return (
    <div className="entity-table">
      {title && <h2 className="entity-table__title">{title}</h2>}
      <div className="entity-table__content">
        {children(tableData)}
      </div>
    </div>
  );
};

export default EntityTable;