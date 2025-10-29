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

export const base64ToBytes = (base64String: string): ArrayBuffer => {
  // converts a base64 encoded string to Uint8Array so it can be parsed by
  // ord-schema js wrapper
  const binaryString = window.atob(base64String);
  const length = binaryString.length;
  const arrayBuffer = new ArrayBuffer(length);
  const view = new Uint8Array(arrayBuffer);
  
  for (let i = 0; i < length; i++) {
    view[i] = binaryString.charCodeAt(i);
  }
  
  return arrayBuffer;
};

export default base64ToBytes;