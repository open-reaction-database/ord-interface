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
import { ActionIcon, Flex, Text, TextInput, Tooltip } from '@mantine/core';
import { IconPlus, IconTrash } from '@tabler/icons-react';
import classes from './SearchItemList.module.scss';

interface SearchItemListProps {
  title: string;
  itemList: string[];
  onUpdateItemList: (newList: string[]) => void;
}

/** Editable list of free-text filter values (reaction IDs, SMARTS, DOIs, …). */
const SearchItemList: React.FC<SearchItemListProps> = ({
  title,
  itemList,
  onUpdateItemList,
}) => {
  const [itemToAdd, setItemToAdd] = useState('');

  const addItem = () => {
    const trimmed = itemToAdd.trim();
    if (!trimmed) return;
    onUpdateItemList([...itemList, trimmed]);
    setItemToAdd('');
  };

  const deleteItem = (idx: number) => {
    onUpdateItemList(itemList.filter((_, index) => index !== idx));
  };

  return (
    <div className={classes.itemList}>
      <Text className={classes.title}>{title}</Text>
      {itemList.map((item, idx) => (
        <Flex
          key={`${item}-${idx}`}
          align="center"
          justify="space-between"
          gap="xs"
          className={classes.row}
        >
          <Tooltip label={item}>
            <Text
              className={classes.rowText}
              truncate
            >
              {item}
            </Text>
          </Tooltip>
          <ActionIcon
            variant="transparent"
            color="secondary.1"
            onClick={() => deleteItem(idx)}
            aria-label={`Remove ${item}`}
          >
            <IconTrash size={16} />
          </ActionIcon>
        </Flex>
      ))}
      <Flex
        gap="xs"
        align="center"
      >
        <TextInput
          className={classes.input}
          size="xs"
          value={itemToAdd}
          onChange={event => setItemToAdd(event.currentTarget.value)}
          onKeyDown={event => {
            if (event.key === 'Enter') {
              event.preventDefault();
              addItem();
            }
          }}
          aria-label={`Add ${title}`}
        />
        <ActionIcon
          variant="default"
          size={30}
          onClick={addItem}
          disabled={!itemToAdd.trim()}
          aria-label={`Add to ${title}`}
        >
          <IconPlus size={16} />
        </ActionIcon>
      </Flex>
    </div>
  );
};

export default SearchItemList;
