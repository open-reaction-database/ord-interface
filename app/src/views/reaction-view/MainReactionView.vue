<script>
import ordSchema from "ord-schema"

export default {
  data() {
    return {
      reaction: {},
      reactionSummary: null,
      loading: true
    }
  },
  computed: {
    reactionId() {
      return this.$route.params.reactionId
    },
  },
  methods: {
    getReactionData () {
      return new Promise(resolve => {
        const xhr = new XMLHttpRequest();
        xhr.open("GET", `/api/getReaction/${this.reactionId}`)
        xhr.responseType = "arraybuffer";
        xhr.onload = () => {
          // if response is good, deserialize reaction and return object from protobuff
          let reaction = null
          if (xhr.response !== null) {
            const bytes = new Uint8Array(xhr.response);
            reaction = ordSchema.reaction_pb.Reaction.deserializeBinary(bytes).toObject();
          }
          resolve(reaction);
        }
        xhr.send()
      })
    },
    getReactionSummary () {
      fetch(`/api/render/${this.reactionId}?compact=false`)
        .then(response => response.json())
        .then(responseData => {
          this.reactionSummary = responseData
        })
    },
  },
  async mounted() {
    this.reaction = await this.getReactionData()
    this.getReactionSummary()
    this.loading = false
  }
}
</script>

<template lang="pug">
.main-reaction-view
  .summary(
    v-if='reactionSummary'
    v-html='reactionSummary'
  )
</template>

<style lang="sass" scoped>

</style>