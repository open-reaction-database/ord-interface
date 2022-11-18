<script>
import EntityTable from '@/components/EntityTable'
import LoadingSpinner from '@/components/LoadingSpinner'

export default {
  props: {
    searchResults: Array,
  },
  components: {
    EntityTable,
    LoadingSpinner,
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
        v-if='row.reactionTable'
      )
      LoadingSpinner(v-else)

</template>

<style lang="sass" scoped>
.search-results-main
  .row
    background-color: white
    border-radius: 0.25rem
    padding: 1rem
    margin-bottom: 1rem
</style>