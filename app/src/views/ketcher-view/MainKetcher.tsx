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

import React, { useEffect, useState } from 'react';
import './MainKetcher.scss';

/**
 * Loads the standalone Ketcher bundle (extracted under app/src/ketcher/ per the
 * project README) into the document. The bundle's main.<hash>.js entry name
 * changes between releases, so resolve it at runtime instead of pinning the
 * hash in the import.
 */
const KETCHER_BASE = '/src/ketcher';

const MainKetcher: React.FC = () => {
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;

    const loadKetcher = async () => {
      try {
        const manifestRes = await fetch(`${KETCHER_BASE}/asset-manifest.json`);
        if (!manifestRes.ok) throw new Error(`asset-manifest.json: HTTP ${manifestRes.status}`);
        const manifest = (await manifestRes.json()) as { files?: Record<string, string> };
        const mainJs = manifest.files?.['main.js'];
        if (!mainJs) throw new Error('asset-manifest.json missing files["main.js"]');

        if (cancelled) return;

        const script = document.createElement('script');
        script.src = `${KETCHER_BASE}/${mainJs.replace(/^\.?\//, '')}`;
        script.async = true;
        script.onerror = () => setError(`Failed to load Ketcher bundle at ${script.src}`);
        document.body.appendChild(script);
      } catch (err) {
        if (!cancelled) {
          console.error('Failed to load Ketcher:', err);
          setError(err instanceof Error ? err.message : 'Failed to load Ketcher');
        }
      }
    };

    loadKetcher();
    return () => {
      cancelled = true;
    };
  }, []);

  return (
    <div id="ketcher">
      <noscript>You need to enable JavaScript to run this app.</noscript>
      {error && <div className="ketcher-load-error">{error}</div>}
      {/* Ketcher's standalone bundle looks up #root to mount itself.
          This route is loaded inside an <iframe> from ModalKetcher, so taking
          over the iframe's React root is intentional. */}
      <div id="root"></div>
    </div>
  );
};

export default MainKetcher;
