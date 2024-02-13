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

import {createRouter, createWebHistory} from 'vue-router'

import About from '../views/About.vue'
import Home from '../views/Home.vue'

const routes =
    [
      {path : '/', name : 'home', component : Home},
      {path : '/about', name : 'about', component : About},
      {
        path : '/browse',
        name : 'browse',
        component : () => import('../views/browse/MainBrowse.vue')
      },
      {
        path : '/browse/set',
        name : 'selected-set',
        component : () =>
            import('../views/browse/selected-set/MainSelectedSet.vue')
      },
      {
        path : '/search',
        name : 'search',
        component : () => import('../views/search/MainSearch.vue')
      },
      {
        path : '/id/:reactionId',
        name : 'reaction-view',
        component : () => import('../views/reaction-view/MainReactionView.vue')
      },
      {
        path : '/ketcher',
        name : 'ketcher',
        component : () => import('../views/viewKetcher/MainKetcher.vue')
      },
    ]

    const router =
        createRouter({history : createWebHistory(process.env.BASE_URL), routes})

export default router
