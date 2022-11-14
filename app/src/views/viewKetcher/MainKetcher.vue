<template lang="pug">
#ketcher
  noscript You need to enable JavaScript to run this app.
  #root
</template>

<script>

export default {
  mounted() {
    // import ketcher js in mounted so the dom is ready
    import("@/views/ketcher/static/js/main.027562ee.js")
  }
}
</script>

<style lang="sass" scoped>
@import '@/views/ketcher/static/css/main.8e693d51.css'
</style>