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

import React, { useState, useEffect, useCallback } from 'react';
import LoadingSpinner from './LoadingSpinner';
import './ModalKetcher.scss';

interface ModalKetcherProps {
  smiles: string;
  onUpdateSmiles: (smiles: string) => void;
  onCloseModal: () => void;
}

const ModalKetcher: React.FC<ModalKetcherProps> = ({ smiles, onUpdateSmiles, onCloseModal }) => {
  const [contWin, setContWin] = useState<Window | null>(null);
  const [mutatedSmiles, setMutatedSmiles] = useState<string>('');
  const [loading, setLoading] = useState<boolean>(true);

  const getKetcher = useCallback(() => {
    const getKetcherInterval = setInterval(() => {
      // attempt to get Ketcher from the iframe once it's loaded
      const iframe = document.getElementById('ketcher-iframe') as HTMLIFrameElement;
      if (!contWin && iframe?.contentWindow) {
        // Check if ketcher is available in the iframe
        try {
          if ((iframe.contentWindow as any).ketcher) {
            // assigning contentWindow.ketcher doesn't seem to persist as expected
            // so we just assign the contentWindow and reference .ketcher
            setContWin(iframe.contentWindow);
            clearInterval(getKetcherInterval);
            drawSmiles();
            setLoading(false);
          }
        } catch (error) {
          // Cross-origin or other access issues - continue trying
          console.log('Waiting for Ketcher to load...');
        }
      }
    }, 1000);

    // Cleanup interval after 30 seconds to avoid infinite polling
    setTimeout(() => {
      clearInterval(getKetcherInterval);
      if (loading) {
        setLoading(false);
        console.warn('Ketcher failed to load within 30 seconds');
      }
    }, 30000);
  }, [contWin, loading]);

  const drawSmiles = useCallback(async () => {
    if (mutatedSmiles && contWin) {
      try {
        // get molblock if we already have SMILES
        const response = await fetch(`/api/molfile?smiles=${encodeURIComponent(mutatedSmiles)}`);
        if (response.ok) {
          const responseData = await response.json();
          (contWin as any).ketcher.setMolecule(responseData);
        } else {
          console.warn('Failed to get molfile for SMILES:', mutatedSmiles);
        }
      } catch (error) {
        console.error('Error loading molecule into Ketcher:', error);
      }
    }
  }, [mutatedSmiles, contWin]);

  const saveSmiles = useCallback(async () => {
    if (contWin) {
      try {
        const newSmiles = await (contWin as any).ketcher.getSmiles();
        setMutatedSmiles(newSmiles);
        onUpdateSmiles(newSmiles);
        onCloseModal();
      } catch (error) {
        console.error('Error getting SMILES from Ketcher:', error);
      }
    }
  }, [contWin, onUpdateSmiles, onCloseModal]);

  const handleBackgroundClick = (e: React.MouseEvent) => {
    // Close modal if clicking on background (not on modal content)
    if (e.target === e.currentTarget) {
      onCloseModal();
    }
  };

  const handleKeyDown = useCallback((e: KeyboardEvent) => {
    // Close modal on Escape key
    if (e.key === 'Escape') {
      onCloseModal();
    }
  }, [onCloseModal]);

  useEffect(() => {
    setMutatedSmiles(smiles);
  }, [smiles]);

  useEffect(() => {
    getKetcher();
  }, [getKetcher]);

  useEffect(() => {
    // Add keyboard event listener
    document.addEventListener('keydown', handleKeyDown);
    return () => {
      document.removeEventListener('keydown', handleKeyDown);
    };
  }, [handleKeyDown]);

  // Re-trigger drawSmiles when contWin becomes available
  useEffect(() => {
    if (contWin && mutatedSmiles) {
      drawSmiles();
    }
  }, [contWin, drawSmiles, mutatedSmiles]);

  return (
    <div className="background" onClick={handleBackgroundClick}>
      <div id="ketcher_modal" className="modal">
        <div className="modal-content">
          <div className="modal-body">
            <iframe
              id="ketcher-iframe"
              src="/ketcher"
              title="Ketcher Molecular Editor"
            />
          </div>
          <div className="modal-footer">
            <button type="button" onClick={onCloseModal}>
              Cancel
            </button>
            <button type="button" onClick={saveSmiles} disabled={loading}>
              Save
            </button>
          </div>
          {loading && (
            <div className="modal-loading">
              <LoadingSpinner />
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default ModalKetcher;