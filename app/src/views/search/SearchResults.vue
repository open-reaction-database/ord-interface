<script>
import EntityTable from '@/components/EntityTable'

export default {
  props: {
    searchResults: Array,
  },
  components: {
    EntityTable
  }
}
</script>

<template lang="pug">
.search-results-main
  EntityTable(
    :tableData='searchResults'
    title="",
    v-slot='{ entities }'
    v-if='searchResults.length'
    :displaySearch='false'
  ) 
    .item(
      v-for='row in entities'
    )
      .id {{row.reaction_id}}

</template>

<style lang="sass" scoped>

</style>