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
      showOptions: true,
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
  .search-options-container(:class='showOptions ? "slide-out" : "hidden"')
    .options-holder
      SearchOptions(
        @searchOptions='updateSearchOptions'
      )
    .slide-out-tab(@click='showOptions=!showOptions')
      .line
      .line
      .line
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
    @media (max-width: 1000px)
      grid-template-columns: 1fr
      .search-options-container
        position: fixed
        height: 100vh
        width: 90%
        transition: 0.5s
        top: 0
        &.hidden
          left: -81%
        &.slide-out
          left: 0
        .options-holder
          padding: 1rem
          background-color: white
          width: 90%
          height: 100%
          overflow-y: auto
          box-shadow: 0 0 5px #a0a0a0
          box-sizing: border-box
        .slide-out-tab
          background-color: white
          position: absolute
          left: 89.5%
          width: 7.5%
          top: 4rem
          box-shadow: 6px 3px 5px #ccc
          padding: 0.5rem 0
          border-bottom-right-radius: 0.25rem
          .line
            background-color: black
            height: 2px
            width: 60%
            margin: 0.5rem auto

</style>