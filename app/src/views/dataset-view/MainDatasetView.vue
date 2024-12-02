<!--
 Copyright 2023 Open Reaction Database Project Authors

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
-->

<script>
import SearchOptions from '../search/SearchOptions'
import SearchResults from './SearchResults'
import ChartView from './ChartView'
import LoadingSpinner from '@/components/LoadingSpinner'
import reaction_pb from "ord-schema"
import base64ToBytes from "@/utils/base64"

export default {
  components: {
    ChartView,
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
      datasetId: "",
      datasetLoading: true,
      datasetData: [],
      isCollapsed: true
    }
  },
  methods: {
    expandOrShrink (e) {
      this.isCollapsed = !this.isCollapsed;
    },
    async getSearchResults() {
      this.loading = true
      // get raw url query string
      this.datasetId =  window.location.pathname.split("/")[2]
      try {
        // If this is the first time we are attempting the search, set the search task ID.
        if (this.searchTaskId == null) {
          // Submit a search task to the server. This returns a GUID task ID.
          const taskres = await fetch(`/api/submit_query?dataset_id=${this.datasetId}&limit=100`, {method: "GET"})
          this.searchTaskId = (await taskres.json());
        }
        // Check the status of the search task.
        const queryRes = await fetch(`/api/fetch_query_result?task_id=${this.searchTaskId}`, {method: "GET"})
        .then(
          (res) => {
            this.searchLoadStatus = res;
            // If one of these codes, search is finished. Return the results or lack of results.
            if (res?.status == 200) {
              this.searchTaskId = null;
              clearInterval(this.searchPollingInterval);
              return res.json();
            }
            else if (res?.status == 400) {
              let taskId = this.searchTaskId;
              this.searchTaskId = null;
              clearInterval(this.searchPollingInterval);
              throw new Error("Error - Search task ID " + taskId + " does not exist");
            }
            else if (res?.status == 500) {
              let taskId = this.searchTaskId;
              this.searchTaskId = null;
              clearInterval(this.searchPollingInterval);
              throw new Error("Error - Search task ID " + taskId + " failed due to server error");
            }
          }
        )
        .then(
          (searchResultsData) => {
            // Deserialize the results.
            this.searchResults = searchResultsData;
            // unpack protobuff for each reaction in results
            this.searchResults.forEach((reaction) => {
                const bytes = base64ToBytes(reaction.proto)
                reaction.data = reaction_pb.Reaction.deserializeBinary(bytes).toObject();
              })
            this.loading = false
          }
        )
      } catch (e) {
        console.log(e)
        this.searchResults = []
        this.loading = false
      }
    },

    // Set Dataset IDs here - 1 or 2
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
      //yield and conversion add if not max values, otherwise remove from query
      if (options.reaction.min_yield !== 0 || options.reaction.max_yield !== 100) {
        this.searchParams["min_yield"] = options.reaction.min_yield
        this.searchParams["max_yield"] = options.reaction.max_yield
      }
      else {
        delete this.searchParams["min_yield"]
        delete this.searchParams["max_yield"]
      }
      if (options.reaction.min_conversion !== 0 || options.reaction.max_conversion !== 100) {
        this.searchParams["min_conversion"] = options.reaction.min_conversion
        this.searchParams["max_conversion"] = options.reaction.max_conversion
      }
      else {
        delete this.searchParams["min_conversion"]
        delete this.searchParams["max_conversion"]
      }

      // general options
      this.searchParams["limit"] = options.general.limit || 100

      // navigate to search page with new params
      this.$router.push({ name: 'search', query: this.searchParams})
    },
  },
  async mounted() {
    this.datasetId =  window.location.pathname.split("/")[2]
    fetch("/api/datasets", {method: "GET"})
      .then(response => response.json())
      .then(data => {
        this.datasetData = data.filter((dataset) => dataset.dataset_id == this.datasetId)[0]
        this.datasetLoading = false
    })
    // Fetch results. If server returns a 102, set up a poll to keep checking back until we have results.
    await this.getSearchResults().then(() =>{
      if (this.searchLoadStatus?.status == 202 && this.searchPollingInterval == null) {
        this.searchPollingInterval = setInterval(() => {this.getSearchResults(), 1000});
        setTimeout(() => {clearInterval(this.searchPollingInterval), 120000});
      }
    })
  },
}
</script>

<template lang="pug">
#dataset-main
  h1 Dataset View
  #charts(:style='this.isCollapsed ? "display: grid" : "display: flex; flex-direction: column"')
    #chartsection(:style='this.isCollapsed ? "width: 40%" : "width: 100%"')
      #chartsectioncharts(:style='this.isCollapsed ? "display: block" : "display: flex"')
        ChartView(
          uniqueId='reactantsFrequency'
          title='Frequency of Reactants'
          apiCall='input_stats'
          role='reactant'
          :isCollapsed='this.isCollapsed'
        )
        ChartView(
          uniqueId='productsFrequency'
          title='Frequency of Products'
          apiCall='product_stats'
          role='product'
          :isCollapsed='this.isCollapsed'
        )
      #expand
        button(@click='expandOrShrink')
          i.material-icons(:style='"margin-top: 15%"' :title='this.isCollapsed ? "Expand" : "Collapse"') {{ this.isCollapsed ? "keyboard_double_arrow_right" : "keyboard_double_arrow_left" }}
    #datasection
      .h4 Dataset Metadata
      table
        tr
          td Dataset ID:
          td {{datasetData.dataset_id ?? '(no id)'}}
        tr
          td Dataset Name:
          td {{datasetData.name ?? '(no name)'}}  
        tr
          td Dataset Description:
          td {{datasetData.description ?? '(no description)'}}
        tr
          td Number of Reactions in Dataset:
          td {{datasetData.num_reactions}}
      .search-results
        SearchResults(
          :searchResults='searchResults'
          :isOverflow='datasetData.num_reactions > searchResults?.length'
          v-if='!loading && searchResults?.length'
        )
        .no-results(v-else-if='!loading && !searchResults?.length')
          .title This dataset contains no reactions.
        .loading(v-else)
          LoadingSpinner
</template>

<style lang="sass" scoped>
@import '@/styles/vars.sass'
#charts
  display: grid
  grid-template-columns: auto 1fr
  column-gap: 1rem
#chartsection
  position: sticky
  width: 100%
  display: flex
  flex-direction: row
#dataset-main
  width: 95%
  margin: 1rem 2.5%
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