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
import reaction_pb from "ord-schema"
import hexToUint from "@/utils/hexToUint"
import ReactionCard from '@/components/ReactionCard'
import DownloadResults from '@/components/DownloadResults'
import CopyButton from '@/components/CopyButton'

export default {
  components: {
    LoadingSpinner,
    ReactionCard,
    DownloadResults,
    CopyButton
  },
  data() {
    return {
      reactions: [],
      loading: true,
      showDownloadResults: false,
    }
  },
  computed: {
    reactionIds() {
      return this.$route.query.reaction_ids || []
    },
    fullUrl() {
      return window.location.href
    }
  },
  methods: {
    async getSelectedReactions() {
      this.loading = true
      try {
        const xhr = new XMLHttpRequest();
        xhr.open("POST", `/api/fetch_reactions`)
        xhr.setRequestHeader("Content-Type", "application/json");
        xhr.onload = () => {
          // if response is good, deserialize reaction and return object from protobuff
          let fetchedReactions = null
          if (xhr.response !== null) {
            fetchedReactions = JSON.parse(xhr.response)
            fetchedReactions.forEach(reaction => {
              const hexString = reaction.proto
              const bytes = hexToUint(hexString)
              reaction.data = reaction_pb.Reaction.deserializeBinary(bytes).toObject();
            })
          }
          this.reactions = fetchedReactions
          this.loading = false
        }
        xhr.send(JSON.stringify(this.reactionIds))
      } catch (e) {
        console.log(e)
        this.loading = false
      }
    },
  },
  mounted() {
    this.getSelectedReactions()
  },
}
</script>

<template lang="pug">
#selected-set-main
  .header
    .title Reaction Set
    .action-button-holder
      CopyButton(
        :textToCopy='fullUrl'
        icon='share'
        buttonText='Shareable Link'
      )
      button(
        :disabled='!reactionIds.length'
        @click='showDownloadResults=true'
      ) Download Reaction Set
  .selected-set(v-if='!loading && reactions.length')
    ReactionCard(
      v-for='reaction in reactions'
      :reaction='reaction'
      :isSelectable='false'
      :isSelected='false'
    )
  .no-results(v-else-if='!loading && !reactions?.length')
    .title There was an issue fetching your selected reactions.
  .loading(v-else)
    LoadingSpinner
  DownloadResults(
    :reactionIds='reactionIds'
    :showDownloadResults='showDownloadResults'
    @hideDownloadResults='showDownloadResults=false'
  )
</template>

<style lang="sass" scoped>
@import '@/styles/vars.sass'
#selected-set-main
  padding: 2rem
  .header
    display: grid
    grid-template-columns: 1fr auto
    margin-bottom: 1rem
    .title
      font-size: 2rem
      font-weight: 700
    .action-button-holder
      display: flex
      align-items: end
      column-gap: 0.5rem
</style>