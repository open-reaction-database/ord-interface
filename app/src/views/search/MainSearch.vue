<script>
import SearchOptions from './SearchOptions'
import SearchResults from './SearchResults'
import LoadingSpinner from '@/components/LoadingSpinner'
import reaction_pb from "ord-schema"
import hexToUint from "@/utils/hexToUint"

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
      showOptions: false,
    }
  },
  methods: {
    async getSearchResults() {
      this.loading = true
      // get raw url query string
      this.urlQuery =  window.location.search
      try {
        const res = await fetch(`/api/query${this.urlQuery}`, {method: "GET"})
        this.searchResults = await res.json()
        // unpack protobuff for each reaction in results
        this.searchResults.forEach((reaction) => {
          const bytes = hexToUint(reaction.proto)
          reaction.data = reaction_pb.Reaction.deserializeBinary(bytes).toObject();
        })
        this.loading = false
      } catch (e) {
        console.log(e)
        this.searchResults = []
        this.loading = false
      }
    },
    updateSearchOptions(options) {
      // reagent options
      if (options.reagent.reagents.length) {
        this.searchParams["component"] = []
        options.reagent.reagents.forEach(reagent => {
          this.searchParams["component"].push(`${reagent.smileSmart};${reagent.source};${reagent.matchMode}`)
        })

        this.searchParams["use_stereochemistry"] = options.reagent.useStereochemistry
        this.searchParams["similarity"] = options.reagent.similarityThreshold
      } else {
        this.searchParams["component"] = []
      }

      // dataset options
      if (options.dataset.datasetIds.length)
        this.searchParams["dataset_ids"] = options.dataset.datasetIds.join(",")
      else
        delete this.searchParams["dataset_ids"]
      if (options.dataset.DOIs.length)
        this.searchParams["dois"] = options.dataset.DOIs.join(",")
      else
        delete this.searchParams["dois"]

      // reaction options
      if (options.reaction.reactionIds.length)
        this.searchParams["reaction_ids"] = options.reaction.reactionIds.join(",")
      else
        delete this.searchParams["reaction_ids"]
      if (options.reaction.reactionSmarts.length)
        this.searchParams["reaction_smarts"] = options.reaction.reactionSmarts.join(",")
      else
        delete this.searchParams["reaction_smarts"]

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
    .title Filters & Options
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
      v-if='!loading && searchResults?.length'
    )
    .no-results(v-else-if='!loading && !searchResults?.length')
      .title No results. Adjust the filters and options and search again.
    .loading(v-else)
      LoadingSpinner

</template>

<style lang="sass" scoped>
@import '@/styles/vars.sass'
#search-main
  width: 95%
  margin: 1rem 2.5%
  display: grid
  grid-template-columns: auto 1fr
  column-gap: 1rem
  min-width: 800px
  .search-options-container
    .title
      font-size: 2rem
      font-weight: 700
      margin-bottom: 0.85rem
    .options-holder
      position: -webkit-sticky
      position: sticky
      top: 1rem
      background-color: white
      padding: 1rem
      box-sizing: border-box
      border-radius: 0.25rem
      overflow-y: auto
  .no-results
    margin-top: 1rem
    text-align: center
  @media (max-width: 1000px)
    grid-template-columns: 1fr
    .title
      display: none
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
        width: 90%
        height: 100%
        overflow-y: auto
        box-shadow: 0 0 5px $darkgrey
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