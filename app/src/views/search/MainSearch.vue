<script>
import SearchOptions from './SearchOptions'
import SearchResults from './SearchResults'

export default {
  components: {
    SearchOptions,
    SearchResults
  },
  data() {
    return {
      searchResults: [],
      queryParams: ""
    }
  },
  methods: {
    updateSearchResults() {
      fetch(`/api/query?${this.queryParams}`, {method: "GET"})
        .then(response => response.json())
        .then(data => {
          console.log('data',data)
          this.searchResults = data
        })
    },
  },
  mounted() {
    // set default query paramaters
    this.queryParams = `dataset_ids=${this.$route.query.datasetId}&limit=100`
    // fetch initial query
    this.updateSearchResults()
  },
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
  SearchResults(
    :searchResults='searchResults'
  )

</template>

<style lang="sass" scoped>
  #search-main
    width: 90%
    margin: 1rem auto
</style>