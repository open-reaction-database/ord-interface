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

import React, { useState } from 'react';
import './CopyButton.scss';

interface CopyButtonProps {
  textToCopy: string;
  icon?: string;
  buttonText: string;
}

const CopyButton: React.FC<CopyButtonProps> = ({ 
  textToCopy, 
  icon = 'copy', 
  buttonText 
}) => {
  const [copied, setCopied] = useState(false);

  const handleCopy = async () => {
    try {
      await navigator.clipboard.writeText(textToCopy);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    } catch (err) {
      console.error('Failed to copy text:', err);
    }
  };

  return (
    <button className="copy-button" onClick={handleCopy}>
      <span className={`copy-button__icon copy-button__icon--${icon}`}>
        {icon === 'share' ? '🔗' : '📋'}
      </span>
      <span className="copy-button__text">
        {copied ? 'Copied!' : buttonText}
      </span>
    </button>
  );
};

export default CopyButton;