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

import React, { useEffect } from 'react';
import './MainKetcher.scss';

const MainKetcher: React.FC = () => {
  useEffect(() => {
    // Import ketcher js in useEffect so the DOM is ready
    // Using dynamic import to load the Ketcher JavaScript
    const loadKetcher = async () => {
      try {
        await import('../../ketcher/static/js/main.027562ee.js' as any);
      } catch (error) {
        console.error('Failed to load Ketcher:', error);
      }
    };

    loadKetcher();
  }, []);

  return (
    <div id="ketcher">
      <noscript>You need to enable JavaScript to run this app.</noscript>
      <div id="root"></div>
    </div>
  );
};

export default MainKetcher;