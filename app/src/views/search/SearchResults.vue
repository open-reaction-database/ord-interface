<script>
export default {
  props: {
    searchResults: Array,
  }
}
</script>

<template lang="pug">
.search-results-main
  .item(
    v-for='item in searchResults'
  )
    .id {{item.dataset_id}}

</template>

<style lang="sass" scoped>

</style>