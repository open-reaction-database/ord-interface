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

export default {
  data() {
    return {
      loading: true, 
      newDatasetName: null,
    }
  },
  methods: {
    async getDatasets() {
      const res = await fetch('/editor-api/getDatasetsByUser/ef4a7184f6bf4d2da99d8e08b36882c0')
      console.log('res.json',res)
    },
    createDataset() {
      if (!this.newDatasetName)
        return alert("Enter a name for the new dataset.")
      // send new dataset name to api
      const xhr = new XMLHttpRequest();
      xhr.datasetName = this.newDatasetName
      xhr.open('POST', `/editor-api/dataset/${this.newDatasetName}/new`);
      xhr.onload = function () {
        if (xhr.status === 409) {
          alert(`Error: Dataset "${this.datasetName}" already exists`);
        } else {
          // get list of datasets
        }
      }
      xhr.send();
    }
  },
  mounted() {
    this.getDatasets()
  }
}
</script>

<template lang="pug">
.create
  .subtitle Create dataset:
  .file-picker
    .input
      label(for='new-dataset-name') Dataset name:
      input#new-dataset-name(
        type='text'
        v-model='newDatasetName'
      )
  .submit
    button#create-dataset(@click='createDataset') Create Dataset
  .copy Lorem ipsum dolor sit amet, consectetur adipiscing elit. Donec tincidunt eros eget magna venenatis, nec volutpat justo ornare. Praesent mauris enim, dignissim nec posuere ut, mollis at lacus. Phasellus ut ultricies arcu, ultrices efficitur enim. Donec pharetra turpis nulla, id elementum felis fermentum sit amet. In id imperdiet sapien. Quisque vulputate odio eget ex aliquet aliquam. Vestibulum a metus magna.
</template>

<style lang="sass" scoped>
.subtitle
  font-size: 1.5rem
.upload
  padding-top: 1rem
  .file-picker
    padding: 1rem 0
    display: flex
    flex-wrap: wrap
    row-gap: 1rem
    .input
      label
        margin-right: 0.5rem
  .submit
    margin-bottom: 2rem
  .copy
    margin-top: 1rem
</style>