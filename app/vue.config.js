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

const {defineConfig} = require('@vue/cli-service')
const path = require('path')

module.exports = defineConfig({
    transpileDependencies: true,
    chainWebpack: config => {
        config.plugin('copy')
            .tap(entries => {
                entries[0].patterns.push({
                    from: path.resolve(__dirname, 'src/ketcher/templates'),
                    to: path.resolve(__dirname, 'dist/templates'),
                    toType: 'dir',
                    noErrorOnMissing: false,
                    globOptions: {ignore: ['.DS_Store']},
                })
                return entries
            })
    },
    devServer: {
      proxy: process.env.NODE_ENV === 'development' ? {
        "^/api": {
          target: "http://0.0.0.0:5000/client",
          changeOrigin: true
        }
      } : {}
    }
})
