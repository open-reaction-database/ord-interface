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

import React from 'react';
import { useLocation } from 'wouter';
import { Button } from '@mantine/core';
import { IconHome } from '@tabler/icons-react';
import classes from './NotFound.module.scss';

interface NotFoundProps {
  code?: number;
  message?: string;
}

// Error state matching ord-app's NotFoundPage; also reused inline for load errors.
const NotFound: React.FC<NotFoundProps> = ({
  code = 404,
  message = 'The page you are looking for does not exist.',
}) => {
  const [, navigate] = useLocation();

  return (
    <div className={classes.container}>
      <div className={classes.code}>{code}</div>
      <div className={classes.message}>{message}</div>
      <Button
        leftSection={<IconHome size={16} />}
        onClick={() => navigate('/')}
      >
        Go to Home Page
      </Button>
    </div>
  );
};

export default NotFound;
