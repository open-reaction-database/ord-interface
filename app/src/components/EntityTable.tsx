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

import React, { useState, useEffect, useMemo, type ReactNode } from 'react';
import './EntityTable.scss';

interface EntityTableProps {
  tableData: any[];
  title?: string;
  displaySearch?: boolean;
  children: (entities: any[]) => ReactNode;
}

const EntityTable: React.FC<EntityTableProps> = ({ 
  tableData, 
  title = '', 
  displaySearch = true,
  children 
}) => {
  const [entities, setEntities] = useState<any[]>([]);
  const [searchString, setSearchString] = useState('');
  const [currentPage, setCurrentPage] = useState(1);
  const [pagination, setPagination] = useState(10);

  useEffect(() => {
    setEntities(tableData);
  }, [tableData]);

  // Reset to first page when pagination changes
  useEffect(() => {
    setCurrentPage(1);
  }, [pagination]);

  const searchArray = useMemo(() => {
    return searchString.split(" ").filter(term => term.trim() !== "");
  }, [searchString]);

  const filteredEntities = useMemo(() => {
    if (!entities) return [];
    let filtered = entities;

    if (searchString && searchArray.length > 0) {
      filtered = entities.filter(entity => {
        // Create a boolean array for each search term
        const matching = Array(searchArray.length).fill(false);
        
        // Check if each search term matches any property in the entity
        searchArray.forEach((param, pIdx) => {
          Object.keys(entity).forEach(key => {
            let value = entity[key];
            if (value == null) value = "null";
            if (value.toString().toLowerCase().includes(param.toLowerCase())) {
              matching[pIdx] = true;
            }
          });
        });
        
        // All search terms must match
        return !matching.includes(false);
      });
    }

    return filtered;
  }, [entities, searchArray]);

  const pagiBottom = useMemo(() => {
    return (currentPage - 1) * pagination;
  }, [currentPage, pagination]);

  const pagiTop = useMemo(() => {
    return currentPage * pagination;
  }, [currentPage, pagination]);

  const paginatedEntities = useMemo(() => {
    return filteredEntities.slice(pagiBottom, pagiTop);
  }, [filteredEntities, pagiBottom, pagiTop]);

  const lastPage = useMemo(() => {
    if (pagination && filteredEntities) {
      return Math.ceil((filteredEntities.length || 1) / pagination);
    }
    return 1;
  }, [pagination, filteredEntities]);

  const pagiPrev = useMemo(() => {
    return currentPage === 1 ? lastPage : currentPage - 1;
  }, [currentPage, lastPage]);

  const pagiNext = useMemo(() => {
    return currentPage === lastPage ? 1 : currentPage + 1;
  }, [currentPage, lastPage]);

  const handlePageChange = (page: number) => {
    if (page >= 1 && page <= lastPage) {
      setCurrentPage(page);
    }
  };

  if (!entities.length) {
    return null;
  }

  return (
    <div className="table-main">
      <div className="table-container">
        <div className="header">
          <div className="title-holder">
            {title && <div className="title">{title}</div>}
          </div>
        </div>
        
        {displaySearch && (
          <div className="search-area">
            <label htmlFor="search">Search: </label>
            <input
              type="text"
              id="search"
              value={searchString}
              onChange={(e) => setSearchString(e.target.value)}
            />
          </div>
        )}
        
        <div className="content">
          <div className="entities-holder">
            {children(paginatedEntities)}
          </div>
          
          {pagination && entities.length > 0 && (
            <div className="pagination">
              <div className="select">
                Showing{' '}
                <select
                  name="pagination"
                  id="pagination"
                  value={pagination}
                  onChange={(e) => setPagination(Number(e.target.value))}
                >
                  <option value={10}>10</option>
                  <option value={25}>25</option>
                  <option value={50}>50</option>
                  <option value={100}>100</option>
                </select>
                {' '}of {filteredEntities.length} entries.
              </div>
              
              <div className="paginav first" onClick={() => handlePageChange(1)}>
                <img className="chevron" src="/img/arrowL.png" alt="Previous" />
                <img className="chevron" src="/img/arrowL.png" alt="Previous" />
                <span className="word">First</span>
              </div>
              
              <div className="prev paginav" onClick={() => handlePageChange(pagiPrev)}>
                <img className="chevron" src="/img/arrowL.png" alt="Previous" />
                <span className="word">Previous</span>
              </div>
              
              <div
                className={`paginav ${currentPage > 1 ? "" : "no-click"}`}
                onClick={() => currentPage > 1 && handlePageChange(currentPage - 1)}
              >
                <span className="word">{currentPage > 1 ? currentPage - 1 : "..."}</span>
              </div>
              
              <div className="paginav no-click">
                <span className="word">{currentPage}</span>
              </div>
              
              <div
                className={`paginav ${currentPage < lastPage ? "" : "no-click"}`}
                onClick={() => currentPage < lastPage && handlePageChange(currentPage + 1)}
              >
                <span className="word">{currentPage < lastPage ? currentPage + 1 : "..."}</span>
              </div>
              
              <div className="next paginav" onClick={() => handlePageChange(pagiNext)}>
                <span className="word">Next</span>
                <img className="chevron" src="/img/arrowR.png" alt="Next" />
              </div>
              
              <div className="paginav last" onClick={() => handlePageChange(lastPage)}>
                <span className="word">Last</span>
                <img className="chevron" src="/img/arrowR.png" alt="Next" />
                <img className="chevron" src="/img/arrowR.png" alt="Next" />
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default EntityTable;