/**
 * Copyright 2026 Open Reaction Database Project Authors
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

// Ported from ord-app's PageContainer/Breadcrumbs.
import { Breadcrumbs as MantineBreadcrumbs, Flex, Tooltip } from '@mantine/core';
import { Link, useLocation } from 'wouter';
import { IconHome } from '@tabler/icons-react';
import classes from './Breadcrumbs.module.scss';

export interface Breadcrumb {
  title: string;
  path: string;
}

interface BreadcrumbsProps {
  items: Breadcrumb[];
}

export function Breadcrumbs({ items }: Readonly<BreadcrumbsProps>) {
  const [, navigate] = useLocation();

  return (
    <div className={classes.container}>
      <MantineBreadcrumbs
        separator="/"
        separatorMargin={6}
        classNames={{ separator: classes.separator, breadcrumb: classes.breadcrumb }}
      >
        {items.map((breadcrumb, index) => {
          const isActive = index === items.length - 1;
          const children = (
            <>
              {index === 0 && (
                <IconHome
                  className={classes.homeIcon}
                  onClick={() => navigate('/')}
                />
              )}
              <Tooltip label={breadcrumb.title}>
                <span className={classes.text}>{breadcrumb.title}</span>
              </Tooltip>
            </>
          );
          return !isActive ? (
            <Link
              className={classes.link}
              to={breadcrumb.path}
              key={breadcrumb.path}
            >
              {children}
            </Link>
          ) : (
            <Flex key={breadcrumb.path}>{children}</Flex>
          );
        })}
      </MantineBreadcrumbs>
    </div>
  );
}

export default Breadcrumbs;
