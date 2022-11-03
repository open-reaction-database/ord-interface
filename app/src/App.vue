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
#app, #main-container
  min-height: 100vh
#main-container
  display: grid
  grid-template-rows: auto 1fr auto
</style>
