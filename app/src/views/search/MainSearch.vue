<script>
import SearchOptions from './SearchOptions'

export default {
  components: {
    SearchOptions
  },
  data() {
    return {
    }
  }
}
</script>

<template lang="pug">
#search-main
  SearchOptions
  button#go_button
    b Search
    //- |     {% if error %}
    //- #error.mt-3.pb-3 {{ error }}
    //- |     {% else %}
    //- #spacer.mt-3.pb-3
    //- |     {% endif %}
  #results.container.mt-3.mb-3
    table#results_table.table

</template>

<style lang="sass" scoped>
  #search-main
    width: 90%
    margin: 1rem auto
</style>