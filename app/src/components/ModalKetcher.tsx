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
import { Button, Flex, Modal } from '@mantine/core';
import { Editor } from 'ketcher-react';
// eslint-disable-next-line @typescript-eslint/ban-ts-comment
// @ts-ignore
import { StandaloneStructServiceProvider } from 'ketcher-standalone';
import type { Ketcher } from 'ketcher-core';
import 'ketcher-react/dist/index.css';
import classes from './ModalKetcher.module.scss';

interface ModalKetcherProps {
  smiles: string;
  onUpdateSmiles: (smiles: string) => void;
  onCloseModal: () => void;
}

const appWindow = globalThis as unknown as { ketcher: Ketcher | null };

const structServiceProvider = new StandaloneStructServiceProvider();

/**
 * Molecule drawing modal on ketcher-react + ketcher-standalone, matching
 * ord-app's ComponentsKetcherEditor. Import lazily — the standalone struct
 * service is a very large chunk.
 */
const ModalKetcher: React.FC<ModalKetcherProps> = ({
  smiles,
  onUpdateSmiles,
  onCloseModal,
}) => {
  const [ketcherInstance, setKetcherInstance] = useState<Ketcher | null>(null);

  useEffect(() => {
    if (ketcherInstance && smiles) {
      ketcherInstance.setMolecule(smiles);
    }
    // Load only the initial structure; subsequent edits live in the ketcher
    // instance itself.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [ketcherInstance]);

  const handleSave = () => {
    if (!ketcherInstance) return;
    ketcherInstance.getSmiles().then(newSmiles => {
      onUpdateSmiles(newSmiles);
      onCloseModal();
    });
  };

  return (
    <Modal
      opened
      onClose={onCloseModal}
      classNames={{ content: classes.wrapper, body: classes.body }}
      withCloseButton={false}
    >
      <div className={classes.editorWrapper}>
        <Editor
          structServiceProvider={structServiceProvider}
          staticResourcesUrl={import.meta.env.BASE_URL as string}
          onInit={ketcher => {
            setKetcherInstance(ketcher);
            // Ketcher EXPECTS global object to have ketcher variable, otherwise it won't work
            appWindow.ketcher = ketcher;
          }}
          errorHandler={console.error}
        />
      </div>
      <Flex
        className={classes.actions}
        gap="sm"
        justify="flex-end"
      >
        <Button
          variant="default"
          onClick={onCloseModal}
        >
          Cancel
        </Button>
        <Button onClick={handleSave}>Save</Button>
      </Flex>
    </Modal>
  );
};

export default ModalKetcher;
