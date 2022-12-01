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
      searchParams: {},
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
      // reagent options
      if (options?.reagent?.reagents?.length) {
        this.searchParams["component"] = []
        options.reagent.reagents.forEach(reagent => {
          this.searchParams["component"].push(`${encodeURIComponent(reagent.smileSmart)};${reagent.source};${reagent.matchMode}`)
        })
        this.searchParams["use_stereochemistry"] = options.reagent.useStereochemistry
        this.searchParams["similarity"] = options.reagent.similarityThreshold
      }

      // dataset options
      if (options.datasetIds.length)
        this.searchParmas["dataset_ids"] = encodeURIComponent(options.datasetIds.join(","))
      if (options.DOIs.length)
        this.searchParams["dois"] = options.DOIs.join(",")

      // reaction options

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