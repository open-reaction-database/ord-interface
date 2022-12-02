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
    },
    downloadResults() {
      // create .pb download of search results
      const xhr = new XMLHttpRequest();
      xhr.open('POST', 'api/download_results');
      xhr.responseType = "blob";
      xhr.onload = () => {
        if (xhr.status === 200) {
          const url = URL.createObjectURL(xhr.response);
          const link = document.createElement('a');
          link.href = url;
          link.download = "ord_search_results.pb.gz"
          link.click();
          // https://stackoverflow.com/a/56547307.
          setTimeout(() => {
            URL.revokeObjectURL(url);
            link.remove();
          }, 100);
        }
      };
      xhr.setRequestHeader('Content-Type', 'application/json');

      // format request json to expected key/value in api
      const requestJson = this.formattedResults.map(result => {return {"Reaction ID": result.reaction_id}})
      xhr.send(JSON.stringify(requestJson));
    },
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
    title="Search Results",
    v-slot='{ entities }'
    v-if='formattedResults.length'
    :displaySearch='false'
  ) 
    .action-button-holder
      button(
        :disabled='!formattedResults.length'
        @click='downloadResults'
      ) Download Results
    template(
      v-for='row in entities'
    )
      router-link(
        :to='{ name: "reaction-view", params: {reactionId: row.reaction_id}}'
      )
        .row
          .id {{row.reaction_id}}
          .reaction-table(
            v-html='row.reactionTable'
            v-if='row.reactionTable'
          )
          LoadingSpinner(v-else)

</template>

<style lang="sass" scoped>
.search-results-main
  margin-top: 2rem
  .action-button-holder
    margin: -2.5rem 0 1rem // bring button row inline with title without having to pass too much into EntityTable
    display: flex
    justify-content: flex-end
    button
  .row
    background-color: white
    border-radius: 0.25rem
    padding: 1rem
    margin-bottom: 1rem
    transition: 0.25s
    &:hover
      box-shadow: 0 0 5px #a0a0a0
      cursor: pointer
</style>