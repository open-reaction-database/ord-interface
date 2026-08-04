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

// Subset of the Ketcher API surface we actually invoke.
interface KetcherApi {
  setMolecule(molfile: string): void;
  getSmiles(): Promise<string>;
}

interface KetcherWindow extends Window {
  ketcher?: KetcherApi;
}

const ModalKetcher: React.FC<ModalKetcherProps> = ({
  smiles,
  onUpdateSmiles,
  onCloseModal,
}) => {
  const [contWin, setContWin] = useState<KetcherWindow | null>(null);
  const [mutatedSmiles, setMutatedSmiles] = useState<string>('');
  const [loading, setLoading] = useState<boolean>(true);

  const drawSmiles = useCallback(async () => {
    if (mutatedSmiles && contWin?.ketcher) {
      try {
        const response = await fetch(
          `/api/molfile?smiles=${encodeURIComponent(mutatedSmiles)}`,
        );
        if (response.ok) {
          const molfile = (await response.json()) as string;
          contWin.ketcher.setMolecule(molfile);
        } else {
          console.warn('Failed to get molfile for SMILES:', mutatedSmiles);
        }
      } catch (error) {
        console.error('Error loading molecule into Ketcher:', error);
      }
    }
  }, [mutatedSmiles, contWin]);

  const saveSmiles = useCallback(async () => {
    if (contWin?.ketcher) {
      try {
        const newSmiles = await contWin.ketcher.getSmiles();
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

  const handleKeyDown = useCallback(
    (e: KeyboardEvent) => {
      // Close modal on Escape key
      if (e.key === 'Escape') {
        onCloseModal();
      }
    },
    [onCloseModal],
  );

  useEffect(() => {
    setMutatedSmiles(smiles);
  }, [smiles]);

  // Poll the iframe's contentWindow once per second until it exposes
  // `ketcher`. Runs once on mount with no React state in the deps array so
  // success doesn't kick off a fresh 30-second idle interval (the previous
  // useCallback-based version re-fired on every relevant state change).
  useEffect(() => {
    const intervalId = window.setInterval(() => {
      const iframe = document.getElementById(
        'ketcher-iframe',
      ) as HTMLIFrameElement | null;
      if (!iframe?.contentWindow) return;
      try {
        const win = iframe.contentWindow as KetcherWindow;
        if (win.ketcher) {
          window.clearInterval(intervalId);
          window.clearTimeout(timeoutId);
          setContWin(win);
          setLoading(false);
        }
      } catch {
        // Cross-origin or other access issues — keep polling.
      }
    }, 1000);

    const timeoutId = window.setTimeout(() => {
      window.clearInterval(intervalId);
      setLoading(prev => {
        if (prev) console.warn('Ketcher failed to load within 30 seconds');
        return false;
      });
    }, 30000);

    return () => {
      window.clearInterval(intervalId);
      window.clearTimeout(timeoutId);
    };
  }, []);

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
    <div
      className="background"
      onClick={handleBackgroundClick}
    >
      <div
        id="ketcher_modal"
        className="modal"
      >
        <div className="modal-content">
          <div className="modal-body">
            {/* Point at Ketcher's own index.html (shipped in the standalone
                bundle under app/public/ketcher/). All of Ketcher's bundled
                asset paths are relative — sourcing from `/ketcher/index.html`
                lets the iframe's document resolve them against
                `/ketcher/static/...` correctly. Pointing the iframe at our
                React Router route instead (the previous design) loaded our
                shell HTML, which broke Ketcher's relative asset resolution. */}
            <iframe
              id="ketcher-iframe"
              src="/ketcher/index.html"
              title="Ketcher Molecular Editor"
            />
          </div>
          <div className="modal-footer">
            <button
              type="button"
              onClick={onCloseModal}
            >
              Cancel
            </button>
            <button
              type="button"
              onClick={saveSmiles}
              disabled={loading}
            >
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
