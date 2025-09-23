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
import { Link } from 'react-router-dom';
import './HeaderNav.scss';

const HeaderNav: React.FC = () => {
  return (
    <nav className="navbar navbar-expand-lg bg-light header-nav">
      <div className="container header-nav__container">
        <a className="navbar-brand header-nav__brand" href="/">
          <img
            src="https://raw.githubusercontent.com/Open-Reaction-Database/ord-schema/main/logos/logo.svg"
            alt="ORD Logo"
            height="30"
          />
        </a>
        <div id="navbarNav" className="collapse navbar-collapse header-nav__nav">
          <div className="navbar-nav header-nav__nav-list">
            <div className="nav-item header-nav__nav-item">
              <Link className="nav-link header-nav__link" to="/browse">Browse</Link>
              <Link className="nav-link header-nav__link" to="/search">Search</Link>
            </div>
            <div className="nav-item header-nav__nav-item">
              <a className="nav-link header-nav__link" href="https://app.open-reaction-database.org">Contribute</a>
              <a className="nav-link header-nav__link" href="https://docs.open-reaction-database.org">Docs</a>
            </div>
            <div className="nav-item header-nav__nav-item">
              <Link className="nav-link header-nav__link" to="/about">About</Link>
            </div>
          </div>
        </div>
      </div>
    </nav>
  );
};

export default HeaderNav;