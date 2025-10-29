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

import React, { useState, useEffect } from 'react';
import { Link } from 'react-router-dom';
import EntityTable from '../../components/EntityTable';
import LoadingSpinner from '../../components/LoadingSpinner';
import './MainBrowse.scss';

interface Dataset {
  dataset_id: string;
  name: string;
  description?: string;
  num_reactions: number;
}

const MainBrowse: React.FC = () => {
  const [loading, setLoading] = useState(true);
  const [tableData, setTableData] = useState<Dataset[]>([]);

  useEffect(() => {
    fetch("/api/datasets", { method: "GET" })
      .then(response => response.json())
      .then((data: Dataset[]) => {
        setTableData(data);
        setLoading(false);
      })
      .catch(error => {
        console.error('Error fetching datasets:', error);
        setLoading(false);
      });
  }, []);

  const truncateDescription = (description?: string) => {
    if (!description) return '';
    return description.length > 75 ? description.substr(0, 75) + "..." : description;
  };

  if (loading) {
    return (
      <div id="browse-main">
        <div className="loading">
          <LoadingSpinner />
        </div>
      </div>
    );
  }

  if (!tableData.length) {
    return (
      <div id="browse-main">
        <div className="loading">
          <LoadingSpinner />
        </div>
      </div>
    );
  }

  return (
    <div id="browse-main">
      <EntityTable tableData={tableData} title="">
        {(entities: Dataset[]) => (
          <div className="table-container">
            <div className="column label">Dataset ID</div>
            <div className="column label">Name</div>
            <div className="column label">Description</div>
            <div className="column label">Size</div>
            {entities.map((row) => (
              <React.Fragment key={row.dataset_id}>
                <div className="column">
                  <Link to={`/dataset/${row.dataset_id}`}>
                    {row.dataset_id}
                  </Link>
                </div>
                <div className="column">{row.name}</div>
                <div className="column">{truncateDescription(row.description)}</div>
                <div className="column">{row.num_reactions}</div>
              </React.Fragment>
            ))}
          </div>
        )}
      </EntityTable>
    </div>
  );
};

export default MainBrowse;