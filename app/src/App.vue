
<script>
import HeaderNav from "./components/HeaderNav"
import MainFooter from "./components/MainFooter"

export default {
  components: {
    HeaderNav,
    MainFooter,
  },
  computed: {
    noHeaderFooter() {
      const noHeaderRoutes = ['ketcher']
      return noHeaderRoutes.includes(this.$route.name)
    },
  }
}
</script>

<template lang="pug">
#main-container(
  :class='noHeaderFooter ? "full-height" : ""'
)
  HeaderNav(v-if='!noHeaderFooter')
  router-view
  MainFooter(v-if='!noHeaderFooter')
</template>

<style lang="sass">
@import url('https://fonts.googleapis.com/css2?family=Roboto:wght@400;700&display=swap')
@import '@/styles/vars'
body
  font-family: 'Roboto', sans-serif
  margin: 0
#app, #main-container
  min-height: 100vh
#main-container
  display: grid
  grid-template-rows: auto 1fr auto
  background-color: $lightgrey
  &.full-height
    grid-template-rows: 1fr
</style>
