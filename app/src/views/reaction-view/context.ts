/**
 * Copyright 2026 Open Reaction Database Project Authors
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

// Read-only analogue of ord-app's reactionContext: the decoded reaction is
// threaded to the section components through context instead of a Redux store.
import { createContext, useContext } from 'react';
import { ord } from 'ord-schema-protobufjs';

export interface ReactionViewContextValue {
  reaction: ord.Reaction;
}

export const reactionViewContext = createContext<ReactionViewContextValue>({
  reaction: new ord.Reaction(),
});

export function useReaction(): ord.Reaction {
  return useContext(reactionViewContext).reaction;
}
