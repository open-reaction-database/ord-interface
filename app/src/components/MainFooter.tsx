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
import { IconBrandGithub, IconBrandLinkedin } from '@tabler/icons-react';
import classes from './MainFooter.module.scss';

const MainFooter: React.FC = () => {
  const year = new Date().getFullYear();

  return (
    <footer className={classes.footer}>
      <div className={classes.inner}>
        <div className={classes.copyright}>
          © Copyright {year} Open Reaction Database
        </div>
        <div>
          <a href="mailto:help@open-reaction-database.org">
            help@open-reaction-database.org
          </a>
        </div>
        <div className={classes.links}>
          <a
            href="https://docs.open-reaction-database.org"
            target="_blank"
            rel="noopener noreferrer"
          >
            Docs
          </a>
          <a
            href="https://github.com/open-reaction-database/"
            target="_blank"
            rel="noopener noreferrer"
            aria-label="GitHub"
          >
            <IconBrandGithub size={18} />
          </a>
          <a
            href="https://www.linkedin.com/company/open-reaction-database/"
            target="_blank"
            rel="noopener noreferrer"
            aria-label="LinkedIn"
          >
            <IconBrandLinkedin size={18} />
          </a>
        </div>
      </div>
    </footer>
  );
};

export default MainFooter;
