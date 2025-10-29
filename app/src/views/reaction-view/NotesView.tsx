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

import React, { useMemo } from 'react';
import './NotesView.scss';

interface NotesViewProps {
  notes: Record<string, any>;
}

interface NoteToDisplay {
  val: any;
  label: string;
}

const NotesView: React.FC<NotesViewProps> = ({ notes }) => {
  const camelToSpaces = (str: string): string => {
    // convert camel case strings to string with spaces
    // also remove leading "is " where needed
    return str.replace(/([a-z])([A-Z])/g, '$1 $2')
             .replace("is ", "")
             .toLowerCase();
  };

  const notesToDisplay = useMemo((): NoteToDisplay[] => {
    if (!notes) return [];
    
    // grab fields with value and that aren't false for display
    return Object.keys(notes)
            .filter(key => notes[key])
            .map(key => ({
              val: notes[key],
              label: camelToSpaces(key)
            }));
  }, [notes]);

  return (
    <div className="notes-view">
      <div className="details">
        {notesToDisplay.map((note, index) => (
          <React.Fragment key={index}>
            <div className="label">{note.label}</div>
            <div className="value">{note.val}</div>
          </React.Fragment>
        ))}
      </div>
    </div>
  );
};

export default NotesView;