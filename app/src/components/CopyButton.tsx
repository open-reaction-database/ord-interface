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
  buttonText?: string;
}

const CopyButton: React.FC<CopyButtonProps> = ({ 
  textToCopy, 
  icon = 'content_copy', 
  buttonText = '' 
}) => {
  const [displayNotification, setDisplayNotification] = useState(false);
  const [notificationStyle, setNotificationStyle] = useState({
    top: 0,
    left: 0
  });

  const handleCopy = async (event: React.MouseEvent<HTMLButtonElement>) => {
    // Move notification block to mouse cursor position
    const { clientX, clientY } = event;
    setNotificationStyle({
      top: clientY,
      left: clientX
    });

    try {
      await navigator.clipboard.writeText(textToCopy);
      setDisplayNotification(true);

      setTimeout(() => {
        setDisplayNotification(false);
      }, 1500);
    } catch (err) {
      console.error('Failed to copy text:', err);
    }
  };

  return (
    <div className="copy-button-main">
      <button onClick={handleCopy}>
        <i className="material-icons">{icon}</i>
        {buttonText && <div className="copy">{buttonText}</div>}
      </button>
      
      {displayNotification && (
        <div 
          id="copy-notification" 
          className={`fade-enter ${displayNotification ? 'fade-enter-active' : ''}`}
          style={{
            top: `${notificationStyle.top}px`,
            left: `${notificationStyle.left}px`
          }}
        >
          Copied to clipboard!
        </div>
      )}
    </div>
  );
};

export default CopyButton;