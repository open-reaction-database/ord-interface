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

/* eslint-disable react-refresh/only-export-components -- utility module; JSX
   is only used for notification icons, nothing here is a component. */
import type { ReactNode } from 'react';
import { notifications, type NotificationData } from '@mantine/notifications';
import { IconCircleCheck, IconCircleX } from '@tabler/icons-react';

// Const object rather than a TS enum: the build runs with erasableSyntaxOnly.
export const NotificationVariant = {
  SUCCESS: 'success',
  ERROR: 'error',
} as const;

export type NotificationVariant =
  (typeof NotificationVariant)[keyof typeof NotificationVariant];

export interface AppNotification extends Omit<NotificationData, 'variant' | 'icon'> {
  variant: NotificationVariant;
}

const iconByVariant: Record<NotificationVariant, ReactNode> = {
  [NotificationVariant.ERROR]: <IconCircleX color="var(--color-red)" />,
  [NotificationVariant.SUCCESS]: <IconCircleCheck color="var(--color-green)" />,
};

/** Toast helper mirroring ord-app's showNotification wrapper. */
export function showNotification({ variant, ...rest }: AppNotification) {
  notifications.show({
    autoClose: 4000,
    icon: iconByVariant[variant],
    color: 'transparent',
    radius: '8px',
    withBorder: false,
    ...rest,
  });
}
