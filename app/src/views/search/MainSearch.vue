<script>
import SearchOptions from './SearchOptions'
import SearchResults from './SearchResults'
import LoadingSpinner from '@/components/LoadingSpinner'

export default {
  components: {
    SearchOptions,
    SearchResults,
    LoadingSpinner
  },
  data() {
    return {
      searchResults: [],
      queryParams: "",
      loading: true
    }
  },
  methods: {
    updateSearchResults() {
      fetch(`/api/query?${this.queryParams}`, {method: "GET"})
        .then(response => response.json())
        .then(data => {
          this.searchResults = data
          this.loading = false
        })
    },
  },
  mounted() {
    // set default query parameters
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
    v-if='!loading'
  )
  .loading(v-else)
    LoadingSpinner

</template>

<style lang="sass" scoped>
  #search-main
    width: 90%
    margin: 1rem auto
</style>