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
import { BrowserRouter as Router, Routes, Route, useLocation } from 'react-router-dom';
import HeaderNav from './components/HeaderNav';
import MainFooter from './components/MainFooter';
import Home from './views/Home';
import About from './views/About';
import MainBrowse from './views/browse/MainBrowse';
import MainSearch from './views/search/MainSearch';
import './App.css';
import './styles/vars.sass';

const AppContent: React.FC = () => {
  const location = useLocation();
  
  // Routes that should not show header/footer
  const noHeaderFooterRoutes = ['ketcher'];
  const noHeaderFooter = noHeaderFooterRoutes.includes(location.pathname.slice(1));

  return (
    <div id="main-container" className={noHeaderFooter ? "full-height" : ""}>
      {!noHeaderFooter && <HeaderNav />}
      <Routes>
        <Route path="/" element={<Home />} />
        <Route path="/about" element={<About />} />
        <Route path="/browse" element={<MainBrowse />} />
        <Route path="/search" element={<MainSearch />} />
      </Routes>
      {!noHeaderFooter && <MainFooter />}
    </div>
  );
};

function App() {
  return (
    <Router>
      <AppContent />
    </Router>
  );
}

export default App;
