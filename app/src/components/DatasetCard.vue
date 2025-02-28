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
import LoadingSpinner from '@/components/LoadingSpinner'
import CopyButton from '@/components/CopyButton'

export default {
  props: {
    reaction: Object,
    isSelected: Boolean,
    isSelectable: Boolean,
    datasetId: String,
    rank: Number
  },
  components: {
    LoadingSpinner,
    CopyButton
  },
  data() {
    return {
      formattedResults: [],
      selectedReactions: [],
      datasetData: []
    }
  },
  methods: {
    goToDataset(event) {
      console.log(event.target)
      // window.location = "/dataset/" + event.target.dataset_id;
    }
  },
  async mounted() {
    fetch("/api/datasets", {method: "GET"})
        .then(response => response.json())
        .then(data => {
          // Returns the dataset with the highest number of reactions by rank, which is passed in.
          this.datasetData = data.sort((a,b) => parseInt(b.num_reactions) - parseInt(a.num_reactions))[this.rank ?? 0]
    })
  },
}
</script>

<template lang="pug">
.dataset-container(v-on:click='goToDataset($event)', dataset_id='')
  .row
    .info
      .col
        .dataset_name {{datasetData?.name}}
        .dataset_description(v-if='datasetData?.description != ""') Description: {{datasetData?.description}}
        .dataset_num_reactions(v-if='datasetData != []') Number of Reactions: {{datasetData?.num_reactions}}
</template>

<style lang="sass" scoped>
@import '@/styles/vars.sass'
.dataset-container
  height: 100px
  width: 23%
  padding: 10px
  font-size: 10pt
  .row
    border-radius: 0.25rem
    padding: 0.75rem
    margin-bottom: 1rem
    transition: 0.25s
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
      .dataset_name
        font-weight: bold
    &.selected
      border-color: $linkblue
</style>
