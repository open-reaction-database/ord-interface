import {createRouter, createWebHistory} from 'vue-router'
import Home from '../views/Home.vue'
import About from '../views/About.vue'

const routes = [
    {
        path: '/',
        name: 'home',
        component: Home
    },
    {
        path: '/about',
        name: 'about',
        component: About
    },
    {
        path: '/browse',
        name: 'browse',
        component: () => import('../views/browse/MainBrowse.vue')
    },
    {
        path: '/browse/set',
        name: 'selected-set',
        component: () => import('../views/browse/selected-set/MainSelectedSet.vue')
    },
    {
        path: '/search',
        name: 'search',
        component: () => import('../views/search/MainSearch.vue')
    },
    {
        path: '/id/:reactionId',
        name: 'reaction-view',
        component: () => import('../views/reaction-view/MainReactionView.vue')
    },
    {
        path: '/ketcher',
        name: 'ketcher',
        component: () => import('../views/viewKetcher/MainKetcher.vue')
    },
]

const router = createRouter({
    history: createWebHistory(process.env.BASE_URL),
    routes
})

export default router
