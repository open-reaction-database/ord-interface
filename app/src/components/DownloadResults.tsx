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

import React, { useEffect, useState } from 'react';
import { Button, Flex, Select, Text } from '@mantine/core';
import { IconDownload } from '@tabler/icons-react';
import FloatingModal from './FloatingModal';
import { api } from '../utils/api';
import { NotificationVariant, showNotification } from '../utils/showNotification';

interface DownloadResultsProps {
  reactionIds: string[];
  showDownloadResults: boolean;
  onHideDownloadResults: () => void;
}

const DownloadResults: React.FC<DownloadResultsProps> = ({
  reactionIds,
  showDownloadResults,
  onHideDownloadResults,
}) => {
  const [fileType, setFileType] = useState<string>('pb.gz');
  const [downloading, setDownloading] = useState(false);

  useEffect(() => {
    const storedFileType = localStorage.getItem('downloadFileType') || 'pb.gz';
    setFileType(storedFileType);
  }, []);

  const handleFileTypeChange = (newFileType: string) => {
    setFileType(newFileType);
    localStorage.setItem('downloadFileType', newFileType);
  };

  const downloadResults = async () => {
    setDownloading(true);
    try {
      const response = await api.post<Blob>(
        '/download_search_results',
        { reaction_ids: reactionIds },
        { responseType: 'blob' },
      );
      const url = URL.createObjectURL(response.data);
      const link = document.createElement('a');
      link.href = url;
      link.download = 'ord_search_results.pb.gz';
      link.click();
      // https://stackoverflow.com/a/56547307.
      setTimeout(() => {
        URL.revokeObjectURL(url);
        link.remove();
      }, 100);
    } catch (error) {
      console.error('Error downloading search results:', error);
      showNotification({
        variant: NotificationVariant.ERROR,
        message: 'Download failed',
      });
    } finally {
      setDownloading(false);
    }
  };

  if (!showDownloadResults) return null;

  return (
    <FloatingModal
      title="Download results"
      size="md"
      onCloseModal={onHideDownloadResults}
    >
      <Flex
        direction="column"
        gap="md"
      >
        <Text size="sm">Select your desired file type and then click download.</Text>
        <Select
          label="File type"
          value={fileType}
          onChange={value => value && handleFileTypeChange(value)}
          allowDeselect={false}
          data={[
            { value: 'pb.gz', label: 'pb.gz' },
            { value: 'csv', label: 'csv (coming soon)', disabled: true },
            { value: 'pbtxt', label: 'pbtxt (coming soon)', disabled: true },
          ]}
        />
        <Flex justify="flex-end">
          <Button
            leftSection={<IconDownload size={16} />}
            loading={downloading}
            onClick={downloadResults}
          >
            Download {fileType} file
          </Button>
        </Flex>
      </Flex>
    </FloatingModal>
  );
};

export default DownloadResults;
