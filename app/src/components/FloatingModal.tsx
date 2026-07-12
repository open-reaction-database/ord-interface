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

import React, { type ReactNode } from 'react';
import { Modal } from '@mantine/core';

interface FloatingModalProps {
  title: ReactNode;
  children: ReactNode;
  onCloseModal: () => void;
  /** Mantine Modal size, e.g. 'md', 'xl', or a CSS width. */
  size?: string | number;
}

// Thin wrapper over Mantine's Modal: callers mount it conditionally, so it is
// always open while rendered.
const FloatingModal: React.FC<FloatingModalProps> = ({
  title,
  children,
  onCloseModal,
  size = 'xl',
}) => {
  return (
    <Modal
      opened
      centered
      radius="sm"
      size={size}
      title={title}
      onClose={onCloseModal}
    >
      {children}
    </Modal>
  );
};

export default FloatingModal;
