<script>
import LoadingSpinner from '@/components/LoadingSpinner'
import reaction_pb from "ord-schema"
import hexToUint from "@/utils/hexToUint"

export default {
  components: {
    LoadingSpinner
  },
  data() {
    return {
      reactionIds: [],
      reactions: [],
      loading: true,
    }
  },
  computed: {
    reactionIds() {
      return this.$route.query.reaction_ids || []
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
            // const hexString = JSON.parse(xhr.response)[0].proto
            // const bytes = hexToUint(hexString)
            // fetchedReactions = reaction_pb.Reaction.deserializeBinary(bytes).toObject();
          }
          this.reactions = fetchedReactions
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
  .selected-set
    .view(v-for)
  .no-results(v-else-if='!loading && !searchResults?.length')
    .title There was an issue fetching your selected reactions.
  .loading(v-else)
    LoadingSpinner

</template>

<style lang="sass" scoped>
@import '@/styles/vars.sass'

</style>