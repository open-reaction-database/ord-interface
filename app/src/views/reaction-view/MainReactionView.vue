<script>
import ordSchema from "ord-schema"

export default {
  data() {
    return {
      reaction: {},
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
        xhr.open("GET", `/api/id/${this.reactionId}`)
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
  },
  async mounted() {
    this.reaction = await this.getReactionData()
  }
}
</script>

<template lang="pug">
.main-reaction-view

</template>

<style lang="sass" scoped>

</style>