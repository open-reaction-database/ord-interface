<script>
import EntityTable from '@/components/EntityTable'

export default {
  props: {
    searchResults: Array,
  },
  components: {
    EntityTable
  },
  data() {
    return {
      formattedResults: []
    }
  },
  methods: {
    getReactionTables() {
      this.formattedResults.forEach(result => {
        fetch(`/api/render/${result.reaction_id}`)
          .then(response => response.json())
          .then(responseData => {
            result.reactionTable = responseData
          })
      })
    }
  },
  async mounted() {
    this.formattedResults = this.searchResults
    this.getReactionTables()
  },
}
</script>

<template lang="pug">
.search-results-main
  EntityTable(
    :tableData='formattedResults'
    title="",
    v-slot='{ entities }'
    v-if='formattedResults.length'
    :displaySearch='false'
  ) 
    .row(
      v-for='row in entities'
    )
      .id {{row.reaction_id}}
      .reaction-table(
        v-html='row.reactionTable'
      )

</template>

<style lang="sass" scoped>

</style>