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

import type {
  Reaction,
  ReactionConditions,
  ReactionNotes,
  ReactionObservation,
  ReactionOutcome,
  ReactionProvenance,
  ReactionSetup,
  ReactionWorkup,
} from 'ord-schema/proto/reaction_pb';

export type ReactionData = Reaction.AsObject;
export type ReactionConditionsData = ReactionConditions.AsObject;
export type ReactionNotesData = ReactionNotes.AsObject;
export type ReactionObservationData = ReactionObservation.AsObject;
export type ReactionOutcomeData = ReactionOutcome.AsObject;
export type ReactionProvenanceData = ReactionProvenance.AsObject;
export type ReactionSetupData = ReactionSetup.AsObject;
export type ReactionWorkupData = ReactionWorkup.AsObject;

export interface SearchResult {
  reaction_id: string;
  // The backend tags each result with the parent dataset id; ReactionCard
  // renders a dataset link from it.
  dataset_id?: string;
  proto: string;
  data: ReactionData;
}

export interface Dataset {
  dataset_id: string;
  name?: string;
  description?: string;
  num_reactions: number;
}

// Mirrors ord_interface.api.nl_query.NLComponent.
export interface NLComponent {
  identifier: string;
  target: 'INPUT' | 'OUTPUT';
  mode: 'EXACT' | 'SIMILAR' | 'SUBSTRUCTURE' | 'SMARTS';
}

// Mirrors ord_interface.api.nl_query.NLQuery: how the model understood the
// question. Component-level constraints plus optional numeric/flag filters.
export interface NLInterpretation {
  components: NLComponent[];
  min_yield?: number | null;
  max_yield?: number | null;
  min_conversion?: number | null;
  max_conversion?: number | null;
  reaction_smarts?: string | null;
  similarity_threshold?: number | null;
  use_stereochemistry?: boolean | null;
  limit?: number | null;
}

// Mirrors ord_interface.api.nl_query.ResolvedComponent: a component after its
// name has been resolved to a concrete SMILES (shown for transparency).
export interface ResolvedComponent {
  identifier: string;
  smiles: string;
  resolver: string;
  target: 'INPUT' | 'OUTPUT';
  mode: 'EXACT' | 'SIMILAR' | 'SUBSTRUCTURE' | 'SMARTS';
}

// Raw /api/nl_query payload, before result protos are deserialized.
export interface NLQueryResponse {
  query: string;
  interpretation: NLInterpretation;
  resolved_components: ResolvedComponent[];
  results: Omit<SearchResult, 'data'>[];
}
