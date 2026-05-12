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

import path from 'node:path'
import {fileURLToPath} from 'node:url'

import vue from '@vitejs/plugin-vue'
import {defineConfig} from 'vite'
import {viteStaticCopy} from 'vite-plugin-static-copy'

const __dirname = path.dirname(fileURLToPath(import.meta.url))

export default defineConfig({
  plugins: [
    vue(),
    // Mirror the chainWebpack copy step from the old vue-cli config: ketcher's
    // chemical templates ship inside src/ketcher/templates and must end up at
    // /templates in the build output.
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
    // The codebase predates Vite and routinely omits `.vue` (and other) file
    // extensions in imports. webpack's default resolver allowed this; Vite
    // requires it to be opted into explicitly.
    extensions: ['.mjs', '.js', '.mts', '.ts', '.jsx', '.tsx', '.json', '.vue'],
  },
  server: {
    port: 8080,
    proxy: {
      '^/api': {
        target: 'http://0.0.0.0:5000',
        changeOrigin: true,
      },
    },
  },
})
