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
import FloatingModal from '../../components/FloatingModal';
import LoadingSpinner from '../../components/LoadingSpinner';
import './MainBrowse.scss';

interface Dataset {
  dataset_id: string;
  name: string;
  description?: string;
  num_reactions: number;
  submitted_at?: string | null;
}

// Description cell: the text is clamped to a few lines, and an always-present
// button opens a floating "card" modal with the dataset's full details. The
// modal overlay means the table layout never reflows.
const DescriptionCell: React.FC<{
  datasetId: string;
  name: string;
  description?: string;
  submittedAt?: string | null;
  numReactions: number;
}> = ({ datasetId, name, description, submittedAt, numReactions }) => {
  const [expanded, setExpanded] = useState(false);

  return (
    <div className="column description-cell">
      <div className="description-text">{description}</div>
      <button
        type="button"
        className="expand-button"
        aria-label="Show dataset details"
        title="Show dataset details"
        onClick={() => setExpanded(true)}
      >
        <svg
          viewBox="0 0 16 16"
          width="14"
          height="14"
          fill="none"
          stroke="currentColor"
          strokeWidth="1.5"
          strokeLinecap="round"
          strokeLinejoin="round"
          aria-hidden="true"
        >
          <path d="M6 2H2v4M10 2h4v4M6 14H2v-4M10 14h4v-4" />
        </svg>
      </button>
      {expanded && (
        <FloatingModal
          title={
            <>
              {name}
              <div className="modal-dataset-id">{datasetId}</div>
            </>
          }
          className="description-modal"
          onCloseModal={() => setExpanded(false)}
        >
          <div className="modal-meta">
            <span>
              <strong>Submitted:</strong> {submittedAt ?? '—'}
            </span>
            <span>
              <strong>Reactions:</strong> {numReactions.toLocaleString()}
            </span>
          </div>
          <div className="modal-description">{description}</div>
        </FloatingModal>
      )}
    </div>
  );
};

const MainBrowse: React.FC = () => {
  const [loading, setLoading] = useState(true);
  const [tableData, setTableData] = useState<Dataset[]>([]);

  useEffect(() => {
    fetch('/api/datasets', { method: 'GET' })
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
      <EntityTable
        tableData={tableData}
        title=""
      >
        {(entities: Dataset[]) => (
          <div className="table-container">
            <div className="column label">Dataset ID</div>
            <div className="column label">Name</div>
            <div className="column label">Description</div>
            <div className="column label size">Size</div>
            {entities.map(row => (
              <React.Fragment key={row.dataset_id}>
                <div className="column">
                  <Link to={`/dataset/${row.dataset_id}`}>{row.dataset_id}</Link>
                </div>
                <div className="column">{row.name}</div>
                <DescriptionCell
                  datasetId={row.dataset_id}
                  name={row.name}
                  description={row.description}
                  submittedAt={row.submitted_at}
                  numReactions={row.num_reactions}
                />
                <div className="column size">{row.num_reactions.toLocaleString()}</div>
              </React.Fragment>
            ))}
          </div>
        )}
      </EntityTable>
    </div>
  );
};

export default MainBrowse;
