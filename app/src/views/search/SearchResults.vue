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
import EntityTable from '@/components/EntityTable'
import LoadingSpinner from '@/components/LoadingSpinner'
import ReactionCard from '@/components/ReactionCard'
import CopyButton from '@/components/CopyButton'
import DownloadResults from '@/components/DownloadResults'

export default {
  props: {
    searchResults: Array,
  },
  components: {
    EntityTable,
    LoadingSpinner,
    CopyButton,
    ReactionCard,
    DownloadResults,
  },
  data() {
    return {
      formattedResults: [],
      selectedReactions: [],
      showDownloadResults: false,
    }
  },
  computed: {
    fullUrl() {
      return window.location.href
    }
  },
  methods: {
    updateSelectedReactions(event) {
      if (event.target.checked) {
        this.selectedReactions.push(event.target.value)
      } else {
        let idx = this.selectedReactions.indexOf(event.target.value)
        if (idx !== -1) {
          this.selectedReactions.splice(idx, 1)
        }
      }
    },
    goToViewSelected() {
      // store storedSet in vuex so we can retrieve it if user comes back from selected-set
      this.$store.commit('setStoredSet', { query: window.location.search, reactions: this.selectedReactions})
      this.$router.push({ name: 'selected-set', query: {reaction_ids: this.selectedReactions}})
    }
  },
  async mounted() {
    this.formattedResults = this.searchResults
    // if query matches storedSet, set selectedReactions
    if (window.location.search == this.$store.state.storedSet?.query) {
      this.selectedReactions = this.$store.state.storedSet.reactions
    }
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
      CopyButton(
        :textToCopy='fullUrl'
        icon='share'
        buttonText='Shareable Link'
      )
      button(
        :disabled='!formattedResults.length'
        @click='showDownloadResults=true'
      ) Download All Search Results
    ReactionCard(
      v-for='row in entities'
      :reaction='row'
      :isSelectable='true'
      :isSelected='selectedReactions.includes(row.reaction_id)'
      @clickedSelect='updateSelectedReactions'
    )
  transition(
    name='fade'
  )
    .view-selected-container(v-if='selectedReactions.length')
      .view-selected-button(@click='goToViewSelected') View {{selectedReactions.length}} selected reactions
  DownloadResults(
    :reactionIds='formattedResults.map(result => result.reaction_id)'
    :showDownloadResults='showDownloadResults'
    @hideDownloadResults='showDownloadResults=false'
  )
</template>

<style lang="sass" scoped>
@import '@/styles/vars.sass'
@import '@/styles/transition.sass'
.search-results-main
  .action-button-holder
    margin: -2.5rem 0 1rem // bring button row inline with title without having to pass too much into EntityTable
    display: flex
    justify-content: flex-end
    column-gap: 0.5rem
    button
  .reaction-container
    text-decoration: none
    .row
      background-color: white
      border-radius: 0.25rem
      padding: 0.75rem
      margin-bottom: 1rem
      transition: 0.25s
      border: 4px solid white
      &:hover
        box-shadow: 0 0 5px $darkgrey
      .select
        input, label
          cursor: pointer
      .reaction-table
        color: black
        overflow-x: wrap
      .info
        display: grid
        grid-template-columns: repeat(2, 50%)
        row-gap: 0.5rem
        column-gap: 1rem
        margin-top: 1rem
        .col.full
          grid-column: 1/3
          button
            font-size: 1.2rem
        .col
          *
            margin-top: 0.5rem
        .smile
          display: flex
          column-gap: 0.5rem
          align-items: center
          margin-top: 0.25rem
          .value
            width: 100%
            white-space: nowrap
            overflow: hidden
            text-overflow: ellipsis
      &.selected
        border-color: $linkblue
  .view-selected-container
    position: fixed
    bottom: 2rem
    right: 6.5rem
    .view-selected-button
      padding: 0.5rem 1rem
      color: white
      background-color: $linkblue
      border-radius: 0.25rem
      cursor: pointer
      box-shadow: 0 0 5px $darkgrey
      animation: bounce 5s infinite ease-in-out
      @keyframes bounce
        0%, 20%
          transform: translateY(0)
        10%
          transform: translateY(-10px)
  @media (max-width: 1000px)
    margin-top: 2.5rem

</style>