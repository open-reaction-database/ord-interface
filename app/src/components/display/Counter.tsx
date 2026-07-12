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

// Ported from ord-app's Counter (common/components/display/Counter).
import { Badge } from '@mantine/core';
import clsx from 'clsx';
import classes from './Counter.module.scss';

interface CounterProps {
  amount: number | string;
  color?: string;
}

export function Counter({ amount, color = 'grey' }: Readonly<CounterProps>) {
  return (
    <Badge
      className={clsx(classes.counter, { [classes.defaultColor]: color === 'grey' })}
      color={color}
    >
      {amount}
    </Badge>
  );
}

export default Counter;
