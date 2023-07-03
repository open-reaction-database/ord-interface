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
  data () {
    return {

    }
  },
  methods: {
    downloadResults () {
      // create .pb download of search results
      const xhr = new XMLHttpRequest();
      xhr.open('POST', 'api/download_results');
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
  }
}
</script>

<template lang="pug">
.download-results-main
  floating-modal(
      v-if='showDownloadResults'
      title="DownloadResults"
      @closeModal='showDownloadResults=false'
    )
      .options Hello World
</template>

<style lang="sass" scoped>
@import '@/styles/vars.sass'
@import '@/styles/transition.sass'
.download-results-main
</style>