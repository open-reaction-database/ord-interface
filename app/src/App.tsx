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
import { Route, Switch } from 'wouter';
import { MantineProvider } from '@mantine/core';
import { Notifications } from '@mantine/notifications';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import HeaderNav from './components/HeaderNav';
import MainFooter from './components/MainFooter';
import NotFound from './views/NotFound';
import Home from './views/Home';
import About from './views/About';
import MainBrowse from './views/browse/MainBrowse';
import MainSelectedSet from './views/browse/selected-set/MainSelectedSet';
import MainSearch from './views/search/MainSearch';
import MainNLSearch from './views/nl-search/MainNLSearch';
import MainDatasetView from './views/dataset-view/MainDatasetView';
import MainReactionView from './views/reaction-view/MainReactionView';
import { theme } from './styles/theme';
import classes from './App.module.scss';

const queryClient = new QueryClient();

const AppContent: React.FC = () => {
  return (
    <div className={classes.shell}>
      <HeaderNav />
      <main className={classes.main}>
        <div className={classes.content}>
          <Switch>
            <Route
              path="/"
              component={Home}
            />
            <Route
              path="/about"
              component={About}
            />
            <Route
              path="/browse"
              component={MainBrowse}
            />
            <Route
              path="/selected-set"
              component={MainSelectedSet}
            />
            <Route
              path="/search"
              component={MainSearch}
            />
            <Route
              path="/ask"
              component={MainNLSearch}
            />
            <Route
              path="/dataset/:datasetId"
              component={MainDatasetView}
            />
            <Route
              path="/id/:reactionId"
              component={MainReactionView}
            />
            <Route>
              <NotFound />
            </Route>
          </Switch>
        </div>
      </main>
      <MainFooter />
    </div>
  );
};

function App() {
  return (
    <MantineProvider theme={theme}>
      <QueryClientProvider client={queryClient}>
        <Notifications
          position="top-right"
          containerWidth={350}
        />
        <AppContent />
      </QueryClientProvider>
    </MantineProvider>
  );
}

export default App;
