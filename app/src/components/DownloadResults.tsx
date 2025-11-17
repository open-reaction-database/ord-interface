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
import FloatingModal from './FloatingModal';
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
  const [fileType, setFileType] = useState<string>('pb.gz');

  useEffect(() => {
    // Get stored file type from localStorage (similar to Vue's $store)
    const storedFileType = localStorage.getItem('downloadFileType') || 'pb.gz';
    setFileType(storedFileType);
  }, []);

  const handleFileTypeChange = (newFileType: string) => {
    setFileType(newFileType);
    // Store selected file type in localStorage
    localStorage.setItem('downloadFileType', newFileType);
  };

  const downloadResults = () => {
    // Create .pb download of search results
    const xhr = new XMLHttpRequest();
    xhr.open('POST', '/api/download_search_results');
    xhr.responseType = 'blob';
    xhr.onload = () => {
      if (xhr.status === 200) {
        const url = URL.createObjectURL(xhr.response);
        const link = document.createElement('a');
        link.href = url;
        link.download = 'ord_search_results.pb.gz';
        link.click();
        // https://stackoverflow.com/a/56547307.
        setTimeout(() => {
          URL.revokeObjectURL(url);
          link.remove();
        }, 100);
      }
    };
    xhr.setRequestHeader('Content-Type', 'application/json');
    xhr.send(JSON.stringify({ reaction_ids: reactionIds }));
  };

  if (!showDownloadResults) return null;

  return (
    <div className="download-results-main">
      <FloatingModal
        title="Download Results"
        onCloseModal={onHideDownloadResults}
      >
        <div className="download-body">
          <div className="title">Select your desired file type and then click download.</div>
          <div className="options">
            <label htmlFor="file-type-select">File type:</label>
            <select
              id="file-type-select"
              value={fileType}
              onChange={(e) => handleFileTypeChange(e.target.value)}
            >
              <option value="pb.gz">pb.gz</option>
              <option value="csv" disabled>csv (coming soon)</option>
              <option value="pbtxt" disabled>pbtxt (coming soon)</option>
            </select>
          </div>
          <div className="download">
            <button onClick={downloadResults}>Download {fileType} file</button>
          </div>
        </div>
      </FloatingModal>
    </div>
  );
};

export default DownloadResults;