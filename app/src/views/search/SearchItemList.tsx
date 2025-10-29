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

import React, { useState, useEffect } from 'react';
import './SearchItemList.scss';

interface SearchItemListProps {
  title: string;
  itemList: string[];
  onUpdateItemList: (newList: string[]) => void;
}

const SearchItemList: React.FC<SearchItemListProps> = ({ title, itemList, onUpdateItemList }) => {
  const [mutatedList, setMutatedList] = useState<string[]>([]);
  const [itemToAdd, setItemToAdd] = useState<string>('');

  const updateList = () => {
    const newList = [...mutatedList, itemToAdd];
    setMutatedList(newList);
    setItemToAdd('');
    onUpdateItemList(newList);
  };

  const deleteItem = (idx: number) => {
    const newList = mutatedList.filter((_, index) => index !== idx);
    setMutatedList(newList);
    onUpdateItemList(newList);
  };

  const handleKeyPress = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter' && itemToAdd.trim()) {
      updateList();
    }
  };

  useEffect(() => {
    setMutatedList(itemList || []);
  }, [itemList]);

  return (
    <div className="item-list">
      <div className="title">{title}</div>
      <div className="list">
        {mutatedList.map((item, idx) => (
          <React.Fragment key={idx}>
            <div className="text">{item}</div>
            <div className="delete">
              <button onClick={() => deleteItem(idx)}>
                <i className="material-icons">delete</i>
              </button>
            </div>
          </React.Fragment>
        ))}
      </div>
      <div className="input">
        <input
          type="text"
          value={itemToAdd}
          onChange={(e) => setItemToAdd(e.target.value)}
          onKeyPress={handleKeyPress}
        />
        <button onClick={updateList} disabled={!itemToAdd.trim()}>
          <i className="material-icons">add</i>
          <span>Add</span>
        </button>
      </div>
    </div>
  );
};

export default SearchItemList;