import { createApp } from 'vue'
import App from './App.vue'
import router from './router'
import { createStore } from 'vuex'

const store = createStore({
  state() {
    return {
      storedSet: null
    }
  },
  mutations: {
    setStoredSet(state, data) {
      state.storedSet = data
    }
  },
  actions: {
    storeSet({ commit }, data) {
      commit('setStoredSet', data)
    }, 
    retrieveSet({ state }) {
      return state.storedSet
    }
  }
})

const app = createApp(App)
app.use(router)
app.use(store)
app.mount('#app')
