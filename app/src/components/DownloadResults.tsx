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

import React from 'react';
import './DownloadResults.scss';

interface DownloadResultsProps {
  reactionIds: string[];
  showDownloadResults: boolean;
  onHideDownloadResults: () => void;
}

const DownloadResults: React.FC<DownloadResultsProps> = ({ 
  reactionIds, 
  showDownloadResults, 
  onHideDownloadResults 
}) => {
  if (!showDownloadResults) return null;

  const handleDownload = (format: string) => {
    // Placeholder download functionality
    console.log(`Downloading ${reactionIds.length} reactions in ${format} format`);
    onHideDownloadResults();
  };

  return (
    <div className="download-results">
      <div className="download-results__overlay" onClick={onHideDownloadResults} />
      <div className="download-results__modal">
        <div className="download-results__header">
          <h3 className="download-results__title">Download Reactions</h3>
          <button 
            className="download-results__close"
            onClick={onHideDownloadResults}
          >
            ×
          </button>
        </div>
        <div className="download-results__content">
          <p>Choose format to download {reactionIds.length} reaction(s):</p>
          <div className="download-results__buttons">
            <button onClick={() => handleDownload('json')}>JSON</button>
            <button onClick={() => handleDownload('csv')}>CSV</button>
            <button onClick={() => handleDownload('xlsx')}>Excel</button>
          </div>
        </div>
      </div>
    </div>
  );
};

export default DownloadResults;