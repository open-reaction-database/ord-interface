import { createApp } from 'vue'
import App from './App.vue'
import router from './router'
import { createStore } from 'vuex'
import PrimeVue from 'primevue/config'

const store = createStore({
  state() {
    return {
      storedSet: null,
      downloadFileType: null,
    }
  },
  mutations: {
    setStoredSet(state, data) {
      state.storedSet = data
    },
    setDownloadFileType(state, data) {
      state.downloadFileType = data
    }
  },
  actions: {
  }
})

const app = createApp(App)
app.use(router)
app.use(store)
app.use(PrimeVue)
app.mount('#app')
