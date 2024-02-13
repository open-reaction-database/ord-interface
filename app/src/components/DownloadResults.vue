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
import FloatingModal from "@/components/FloatingModal"

export default {
  name: 'DownloadResults',
  props: {
    reactionIds: Array,
    showDownloadResults: Boolean,
  },
  components: {
    "floating-modal": FloatingModal
  },
  watch: {
    fileType(val) {
      // store selected file type in vue session
      this.$store.commit('setDownloadFileType', val)
    }
  },
  data () {
    return {
      fileType: "pb.gz"
    }
  },
  methods: {
    downloadResults () {
      // create .pb download of search results
      const xhr = new XMLHttpRequest();
      xhr.open('POST', '/api/download_results');
      xhr.responseType = "blob";
      xhr.onload = () => {
        if (xhr.status === 200) {
          const url = URL.createObjectURL(xhr.response);
          const link = document.createElement('a');
          link.href = url;
          link.download = "ord_search_results.pb.gz"
          link.click();
          // https://stackoverflow.com/a/56547307.
          setTimeout(() => {
            URL.revokeObjectURL(url);
            link.remove();
          }, 100);
        }
      };
      xhr.setRequestHeader('Content-Type', 'application/json');

      // format request json to expected key/value in api
      const requestJson = this.reactionIds.map(reactionId => {return {"Reaction ID": reactionId}})
      xhr.send(JSON.stringify(requestJson));
    }
  },
  mounted() {
    this.fileType = this.$store.state.downloadFileType || "pb.gz"
  }
}
</script>

<template lang="pug">
.download-results-main
  floating-modal(
      v-if='showDownloadResults'
      title="DownloadResults"
      @closeModal='$emit("hideDownloadResults")'
    )
      .download-body
        .title Select your desired file type and then click download.
        .options 
          label(for='file-type-select') File type:
          select#file-type-select(v-model='fileType')
            option(value="pb.gz") pb.gz
            option(value="csv" disabled) csv (coming soon)
            option(value="pbtxt" disabled) pbtxt (coming soon)
        .download
          button(@click='downloadResults') Download {{fileType}} file
</template>

<style lang="sass" scoped>
@import '@/styles/vars.sass'
@import '@/styles/transition.sass'
.download-body
  height: 100%
  display: grid
  grid-template-rows: 50% auto 1fr
  *
    margin: auto
    width: fit-content
  .title
    font-size: 1.5rem
  .options
    label
      margin-right: 0.5rem
    select
      min-width: 7rem
  .download
    margin-top: 2rem
</style>