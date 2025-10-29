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
import './FloatingModal.scss';

interface FloatingModalProps {
  title: string;
  children: ReactNode;
  onCloseModal: () => void;
}

const FloatingModal: React.FC<FloatingModalProps> = ({
  title,
  children,
  onCloseModal
}) => {
  return (
    <div className="modal-main">
      <div className="modal-background" onClick={onCloseModal} />
      <div className="modal-holder">
        <div className="modal-container">
          <div className="header">
            <div className="title">{title}</div>
            <div className="close" onClick={onCloseModal}>&#10005;</div>
          </div>
          <div className="body">
            {children}
          </div>
        </div>
      </div>
    </div>
  );
};

export default FloatingModal;