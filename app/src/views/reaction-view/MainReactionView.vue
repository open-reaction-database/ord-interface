<script>
import { reaction_pb } from "ord-schema"

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
            reaction = reaction_pb.Reaction.deserializeBinary(bytes).toObject();
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
    getReactionType (id) {
      const identifiers = reaction_pb.ReactionIdentifier.ReactionIdentifierType
      return Object.keys(identifiers).find(key => identifiers[key] == id)
    }
  },
  async mounted() {
    this.reaction = await this.getReactionData()
    this.getReactionSummary()
    this.loading = false
    console.log('schema',reaction_pb)
  }
}
</script>

<template lang="pug">
.main-reaction-view
  .section.summary(v-if='reactionSummary')
    .display(v-html='reactionSummary')
  .section(v-if='reaction?.identifiersList?.length')
    .title Identifiers
    .identifiers
      template(v-for='identifier in reaction.identifiersList')
        .value {{getReactionType(identifier.type)}}
        .value {{identifier.value}}
        .value {{identifier.details}}

</template>

<style lang="sass" scoped>
.section
  width: calc(90vw)
  background-color: white
  border-radius: 0.25rem
  margin: 1rem auto 0
  padding: 1rem
  &.summary
    overflow-x: scroll
</style>