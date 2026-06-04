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
import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import HeaderNav from './components/HeaderNav';
import MainFooter from './components/MainFooter';
import Home from './views/Home';
import About from './views/About';
import MainBrowse from './views/browse/MainBrowse';
import MainSelectedSet from './views/browse/selected-set/MainSelectedSet';
import MainSearch from './views/search/MainSearch';
import MainDatasetView from './views/dataset-view/MainDatasetView';
import MainReactionView from './views/reaction-view/MainReactionView';
import './App.scss';

const queryClient = new QueryClient();

const AppContent: React.FC = () => {
  return (
    <div id="main-container">
      <HeaderNav />
      <Routes>
        <Route
          path="/"
          element={<Home />}
        />
        <Route
          path="/about"
          element={<About />}
        />
        <Route
          path="/browse"
          element={<MainBrowse />}
        />
        <Route
          path="/selected-set"
          element={<MainSelectedSet />}
        />
        <Route
          path="/search"
          element={<MainSearch />}
        />
        <Route
          path="/dataset/:datasetId"
          element={<MainDatasetView />}
        />
        <Route
          path="/id/:reactionId"
          element={<MainReactionView />}
        />
      </Routes>
      <MainFooter />
    </div>
  );
};

function App() {
  return (
    <QueryClientProvider client={queryClient}>
      <Router>
        <AppContent />
      </Router>
    </QueryClientProvider>
  );
}

export default App;
