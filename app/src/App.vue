<template lang="pug">
#main-container
  HeaderNav
  router-view
  MainFooter
</template>

<script>
import HeaderNav from "./components/HeaderNav"
import MainFooter from "./components/MainFooter"

export default {
  components: {
    HeaderNav,
    MainFooter,
  }
}
</script>

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
</style>
