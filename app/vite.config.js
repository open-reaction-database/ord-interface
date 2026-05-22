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

import path from 'node:path';
import { fileURLToPath } from 'node:url';

import vue from '@vitejs/plugin-vue';
import { defineConfig } from 'vite';
import { viteStaticCopy } from 'vite-plugin-static-copy';

const __dirname = path.dirname(fileURLToPath(import.meta.url));

export default defineConfig({
  plugins: [
    vue(),
    // Ketcher's chemical templates ship inside src/ketcher/templates and must
    // end up at /templates in the build output (the Ketcher JS bundle fetches
    // them at runtime). Ketcher v2.5.1 is otherwise self-contained — its main
    // JS bundles all assets except an `overlay.svg` referenced via a relative
    // url() that Vite resolves at build time. If Ketcher is ever updated to a
    // chunked build, the copy targets here will need to grow.
    viteStaticCopy({
      targets: [
        {
          src: 'src/ketcher/templates/*',
          dest: 'templates',
        },
      ],
    }),
  ],
  resolve: {
    alias: {
      '@': path.resolve(__dirname, 'src'),
    },
  },
  server: {
    port: 8080,
    proxy: {
      // 127.0.0.1 (not localhost / 0.0.0.0): on macOS, localhost resolves to
      // ::1 first, and AirPlay Receiver also listens on *:5000. uvicorn binds
      // 127.0.0.1 only, so the explicit IPv4 address dodges the collision.
      '^/api': {
        target: 'http://127.0.0.1:5000',
        changeOrigin: true,
      },
    },
  },
});
