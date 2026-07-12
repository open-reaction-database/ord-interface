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

import { ord } from 'ord-schema-protobufjs';
import { Buffer } from 'buffer';

/**
 * Decodes a base64-encoded Reaction protobuf from the API, matching ord-app's
 * parseReaction (ui/src/store/entities/reactions/reactions.utils.ts).
 */
export function decodeReaction(binpb: string): ord.Reaction {
  return ord.Reaction.decode(Buffer.from(binpb, 'base64'));
}

/** Encodes a Compound message to bytes, e.g. for the /api/compound_svg endpoint. */
export function encodeCompound(compound: ord.ICompound): ArrayBuffer {
  // .finish() returns a view into protobufjs's shared pool buffer; copy to an
  // exact-size buffer because axios transmits the view's whole underlying
  // ArrayBuffer.
  return ord.Compound.encode(compound).finish().slice().buffer;
}

/**
 * Sorted [name, input] entries of Reaction.inputs, ordered by additionOrder.
 * protobufjs represents the inputs map as a plain object, so ordering has to
 * be applied at read time.
 */
export function sortedInputEntries(
  inputs: { [k: string]: ord.IReactionInput } | null | undefined,
): Array<[string, ord.IReactionInput]> {
  return Object.entries(inputs ?? {}).sort(
    (a, b) => (a[1].additionOrder ?? 0) - (b[1].additionOrder ?? 0),
  );
}
