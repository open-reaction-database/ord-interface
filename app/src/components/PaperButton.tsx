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

import React, { type CSSProperties, type ReactNode } from 'react';
import { Button } from '@mantine/core';
import classes from './PaperButton.module.scss';

interface PaperButtonProps {
  title: string;
  description: string;
  icon: ReactNode;
  /** CSS color value for the accent, e.g. 'var(--color-blue)'. */
  color: string;
  /** Renders as a link when set; use onClick otherwise. */
  href?: string;
  onClick?: () => void;
}

/**
 * Action tile with a colored icon chip that inverts to a solid fill on hover.
 * Ported from ord-app's PaperButton so the two frontends share the pattern.
 */
const PaperButton: React.FC<PaperButtonProps> = ({
  title,
  description,
  icon,
  color,
  href,
  onClick,
}) => {
  const style = { '--paper-button-color': color } as CSSProperties;
  const buttonProps = href ? ({ component: 'a', href } as const) : ({} as const);

  return (
    <Button
      {...buttonProps}
      variant="default"
      radius="sm"
      style={style}
      classNames={{
        root: classes.root,
        inner: classes.inner,
        section: classes.section,
        label: classes.label,
      }}
      leftSection={icon}
      onClick={onClick}
    >
      <span>{title}</span>
      <span className={classes.subtitle}>{description}</span>
    </Button>
  );
};

export default PaperButton;
