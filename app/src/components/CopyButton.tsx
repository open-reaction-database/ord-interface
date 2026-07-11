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
import { ActionIcon, Button, Tooltip } from '@mantine/core';
import { IconCopy, IconShare } from '@tabler/icons-react';
import { NotificationVariant, showNotification } from '../utils/showNotification';

interface CopyButtonProps {
  textToCopy: string;
  icon?: 'copy' | 'share';
  buttonText?: string;
}

const CopyButton: React.FC<CopyButtonProps> = ({
  textToCopy,
  icon = 'copy',
  buttonText = '',
}) => {
  const handleCopy = async () => {
    try {
      await navigator.clipboard.writeText(textToCopy);
      showNotification({
        variant: NotificationVariant.SUCCESS,
        message: 'Copied to clipboard',
      });
    } catch (err) {
      console.error('Failed to copy text:', err);
      showNotification({
        variant: NotificationVariant.ERROR,
        message: 'Failed to copy to clipboard',
      });
    }
  };

  const iconNode = icon === 'share' ? <IconShare size={16} /> : <IconCopy size={16} />;

  if (!buttonText) {
    return (
      <Tooltip label="Copy to clipboard">
        <ActionIcon
          variant="transparent"
          color="secondary.1"
          onClick={handleCopy}
          aria-label="Copy to clipboard"
        >
          {iconNode}
        </ActionIcon>
      </Tooltip>
    );
  }

  return (
    <Button
      variant="default"
      radius="sm"
      leftSection={iconNode}
      onClick={handleCopy}
    >
      {buttonText}
    </Button>
  );
};

export default CopyButton;
