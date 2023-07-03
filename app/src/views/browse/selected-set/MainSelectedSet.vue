<script>
import LoadingSpinner from '@/components/LoadingSpinner'
import reaction_pb from "ord-schema"
import hexToUint from "@/utils/hexToUint"
import ReactionCard from '@/components/ReactionCard'

export default {
  components: {
    LoadingSpinner,
    ReactionCard
  },
  data() {
    return {
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
</template>

<style lang="sass" scoped>
@import '@/styles/vars.sass'
#selected-set-main
  padding: 2rem
</style>