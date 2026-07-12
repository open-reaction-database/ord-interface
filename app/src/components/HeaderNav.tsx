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
import { Link, useLocation } from 'wouter';
import { Button } from '@mantine/core';
import { IconExternalLink } from '@tabler/icons-react';
import clsx from 'clsx';
import Wordmark from '../assets/ord-wordmark.svg?react';
import classes from './HeaderNav.module.scss';

const NAV_ITEMS = [
  { label: 'Browse', href: '/browse' },
  { label: 'Search', href: '/search' },
  // The /ask route exists but is intentionally unlisted while the
  // natural-language feature is in development.
  { label: 'About', href: '/about' },
];

const HeaderNav: React.FC = () => {
  const [location] = useLocation();

  return (
    <header className={classes.header}>
      <div className={classes.inner}>
        <Link
          href="/"
          className={classes.brand}
          aria-label="Open Reaction Database home"
        >
          <Wordmark className={classes.wordmark} />
        </Link>

        <nav className={classes.nav}>
          {NAV_ITEMS.map(item => (
            <Link
              key={item.href}
              href={item.href}
              className={clsx(classes.link, {
                [classes.active]: location.startsWith(item.href),
              })}
            >
              {item.label}
            </Link>
          ))}
          <a
            className={classes.link}
            href="https://docs.open-reaction-database.org"
            target="_blank"
            rel="noopener noreferrer"
          >
            Docs
          </a>
        </nav>

        <div className={classes.actions}>
          <Button
            component="a"
            href="https://app.open-reaction-database.org"
            size="compact-sm"
            radius="sm"
            rightSection={<IconExternalLink size={14} />}
          >
            Contribute
          </Button>
        </div>
      </div>
    </header>
  );
};

export default HeaderNav;
