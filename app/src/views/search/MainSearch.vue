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
  watch: {
    '$route.query': {
      handler() {
        this.getSearchResults()
      },
      deep: true,
    }
  },
  data() {
    return {
      searchResults: [],
      queryParams: "",
      searchParams: {},
      loading: true,
      urlQuery: "",
    }
  },
  methods: {
    getSearchResults() {
      this.loading = true
      // get raw url query string
      this.urlQuery =  window.location.search

      fetch(`/api/query${this.urlQuery}`, {method: "GET"})
        .then(response => response.json())
        .then(data => {
          this.searchResults = data
          this.loading = false
        })
    },
    updateSearchOptions(options) {
      // reagent options
      if (options.reagent.reagents.length) {
        this.searchParams["component"] = []
        options.reagent.reagents.forEach(reagent => {
          this.searchParams["component"].push(`${encodeURIComponent(reagent.smileSmart)};${reagent.source};${reagent.matchMode}`)
        })
        this.searchParams["use_stereochemistry"] = options.reagent.useStereochemistry
        this.searchParams["similarity"] = options.reagent.similarityThreshold
      }

      // dataset options
      if (options.dataset.datasetIds.length)
        this.searchParams["dataset_ids"] = encodeURIComponent(options.dataset.datasetIds.join(","))
      if (options.dataset.DOIs.length)
        this.searchParams["dois"] = encodeURIComponent(options.dataset.DOIs.join(","))

      // reaction options
      if (options.reaction.reactionIds.length)
        this.searchParams["reaction_ids"] = encodeURIComponent(options.reaction.reactionIds.join(","))
      if (options.reaction.reactionSmarts.length)
        this.searchParams["reaction_smarts"] = encodeURIComponent(options.reaction.reactionSmarts.join(","))

      // general options
      this.searchParams["limit"] = options.general.limit || 100

      // navigate to search page with new params
      this.$router.push({ name: 'search', query: this.searchParams})
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
  .search-options
    SearchOptions(
      @searchOptions='updateSearchOptions'
    )
  .search-results
    SearchResults(
      :searchResults='searchResults'
      v-if='!loading'
    )
    .loading(v-else)
      LoadingSpinner

</template>

<style lang="sass" scoped>
  #search-main
    width: 95%
    margin: 1rem 2.5%
    display: grid
    grid-template-columns: auto 1fr
    column-gap: 1rem
    min-width: 800px
    .search-options
      position: sticky
      top: 1rem
    .search-results
</style>