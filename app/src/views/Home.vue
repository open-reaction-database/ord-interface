<script>
export default {
  name: 'HomeView',
  components: {
  }
}
</script>

<template lang="pug">
.home
  .home-img
    img( alt="Ord logo" src="https://raw.githubusercontent.com/Open-Reaction-Database/ord-schema/main/logos/logo.svg")
</template>

<style lang="sass" scoped>
.home
  margin-top: 20vh
  .home-img
    width: 90%
    margin: auto
    max-width: 800px
    img
      width: 100%
</style>
