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

import { describe, expect, it } from 'vitest';
import { base64ToBytes } from './base64';

const decode = (base64String: string): number[] =>
  Array.from(new Uint8Array(base64ToBytes(base64String)));

const encode = (bytes: number[]): string => btoa(String.fromCharCode(...bytes));

describe('base64ToBytes', () => {
  it('decodes to a buffer of the right length', () => {
    const buffer = base64ToBytes(btoa('hello'));
    expect(buffer).toBeInstanceOf(ArrayBuffer);
    expect(buffer.byteLength).toBe(5);
  });

  it('decodes ASCII', () => {
    expect(decode(btoa('hello'))).toEqual([104, 101, 108, 108, 111]);
  });

  it('decodes an empty string', () => {
    expect(decode('')).toEqual([]);
  });

  // Serialized protos are arbitrary bytes, so anything that truncated to 7 bits
  // or tripped over a NUL would corrupt every reaction on the page.
  it('round-trips the full byte range', () => {
    const bytes = Array.from({ length: 256 }, (_, i) => i);
    expect(decode(encode(bytes))).toEqual(bytes);
  });

  it('decodes padded and unpadded lengths alike', () => {
    expect(decode(encode([1]))).toEqual([1]);
    expect(decode(encode([1, 2]))).toEqual([1, 2]);
    expect(decode(encode([1, 2, 3]))).toEqual([1, 2, 3]);
  });
});
