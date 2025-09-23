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
import './MainFooter.scss';

const MainFooter: React.FC = () => {
  return (
    <footer className="main-footer">
      <div className="main-footer__container">
        <div className="main-footer__content">
          <p className="main-footer__text">
            © 2024 Open Reaction Database Project
          </p>
          <div className="main-footer__links">
            <a 
              href="https://github.com/Open-Reaction-Database" 
              className="main-footer__link"
              target="_blank" 
              rel="noopener noreferrer"
            >
              GitHub
            </a>
            <a 
              href="https://docs.open-reaction-database.org" 
              className="main-footer__link"
              target="_blank" 
              rel="noopener noreferrer"
            >
              Documentation
            </a>
            <img
              src="https://raw.githubusercontent.com/Open-Reaction-Database/ord-schema/main/logos/logo.svg"
              alt="ORD Logo"
              className="main-footer__logo"
            />
          </div>
        </div>
      </div>
    </footer>
  );
};

export default MainFooter;