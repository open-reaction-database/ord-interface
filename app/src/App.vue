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
</style>
