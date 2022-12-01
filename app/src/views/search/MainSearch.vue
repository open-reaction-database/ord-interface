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
      searchParams: "",
      loading: true
    }
  },
  computed: {
    urlQuery() {
      // get raw url query string
      return window.location.search
    }
  },
  methods: {
    getSearchResults() {
      fetch(`/api/query${this.urlQuery}`, {method: "GET"})
        .then(response => response.json())
        .then(data => {
          this.searchResults = data
          this.loading = false
        })
    },
    updateSearchOptions(options) {
      console.log('options',options)

    },
  },
  mounted() {
    // fetch initial query
    this.getSearchResults()
  },
}
</script>

<template lang="pug">
#search-main
  SearchOptions(
    @searchOptions='updateSearchOptions'
  )
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