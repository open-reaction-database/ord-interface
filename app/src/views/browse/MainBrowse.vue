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
import LoadingSpinner from '../../components/LoadingSpinner.vue'
import MultiRangeSlider from 'multi-range-slider-vue';

export default {
  components: {
    EntityTable,
    LoadingSpinner,
    MultiRangeSlider
  },
  data() {
    return {
      loading: true,
      tableData: []
    }
  },
  mounted() {
    fetch("/api/datasets", {method: "GET"})
      .then(response => response.json())
      .then(data => {
        this.tableData = data
        this.loading = false
      })
  }
}
</script>

<template lang="pug">
#browse-main
  #filters
    label(for='mined') Hide Mined Datasets?
    input#mined(
      type='checkbox'
      checked=true
    )
    label(for='size') Dataset Size (Number of Reactions)
    MultiRangeSlider(
      :miin='0'
      :max='100'
      :ruler='false'
      :label='false'
    )
  EntityTable(
    :tableData='tableData'
    title="",
    v-slot='{ entities }'
    v-if='tableData.length'
  ) 
    .table-container
      .column.label Dataset ID
      .column.label Name
      .column.label Description
      .column.label Size
      template(
        v-for='(row, idx) in entities'
      )
        .column 
          router-link(:to='{ name: "dataset-view", params: {datasetId: row.dataset_id}}') {{row.dataset_id}}
        .column {{row.name}}
        .column {{row.description?.length > 75 ? row.description.substr(0,75)+"..." : row.description}}
        .column {{row.num_reactions}}
  .loading(v-else)
    LoadingSpinner
</template>

<style lang="sass" scoped>
@import '@/styles/table'
#browse-main
  padding: 1rem 0
  .table-container
    grid-template-columns: 1fr 1fr 1fr auto
  .loading
    margin-top: 30vh
#filters
  padding: 10px
  width: 200px
</style>
