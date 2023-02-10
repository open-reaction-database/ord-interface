<script>
import { reaction_pb } from "ord-schema"

export default {
  data() {
    return {
      reaction: {},
      reactionSummary: null,
      reactionBytes: null,
      loading: true,
      inputsIdx: 0,
    }
  },
  computed: {
    reactionId() {
      return this.$route.params.reactionId
    },
    displayInputs() {
      if (!this.reaction) return {}
      let returnArr = {...this.reaction.inputsMap[this.inputsIdx][1]}
      // filter out null/undefined values and arrays
      return Object.fromEntries(Object.entries(returnArr).filter(([_,v]) => v != null && !Array.isArray(v)))
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
            this.reactionBytes = new Uint8Array(xhr.response);
            reaction = reaction_pb.Reaction.deserializeBinary(this.reactionBytes).toObject();
            // sort inputs by addition order
            reaction.inputsMap.sort((a,b) => a[1].additionOrder - b[1].additionOrder)
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
    },
    getCompoundSVG (component) {
      // prepare compound for http request
      const compoundStr = component[1]
                        .componentsList[0]
                        .identifiersList
                        .find((identifier) => identifier.type == 2)
                        .value
      const compound = new reaction_pb.Compound()
      const identifier = compound.addIdentifiers()
      identifier.setValue(compoundStr)
      identifier.setType(reaction_pb.CompoundIdentifier.CompoundIdentifierType.SMILES)
      // send http request
      return new Promise(resolve => {
        const xhr = new XMLHttpRequest()
        xhr.open("POST", "/api/render/compound")
        const binary = compound.serializeBinary()
        xhr.responseType = "json"
        xhr.onload = function () {
          console.log(xhr.response)
          if (xhr.response !== null) {
            console.log("success!")
          }
          resolve(xhr.response)
        }
        xhr.send(binary)
      })
      // return compound
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
  .section(v-if='reaction?.inputsMap?.length')
    .title Inputs
    .tabs
      .tab(
        v-for='(input, idx) in reaction.inputsMap'
        @click='inputsIdx = idx'
        :class='inputsIdx == idx ? "selected" : ""'
      ) {{input[0]}}
    .input
      .title Details
      .details
        template(v-for='key in Object.keys(displayInputs)')
          .label {{key.replaceAll(/(?=[A-Z])/g, ' ')}}
          .value {{displayInputs[key]}}
      .title Components
      .details 
        .svg(
          v-html='getCompoundSVG(reaction.inputsMap[inputsIdx])'
        )
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
  .title
    font-weight: 700
    font-size: 1.5rem
    margin-bottom: 0.5rem
  .tabs
    display: flex
    column-gap: 0.5rem
    margin-bottom: 0.5rem
    .tab
      padding: 0.5rem 1rem
      border-radius: 0.25rem
      border: 1px solid lightgrey
      cursor: pointer
      transition: 0.25s
      &.selected
        background-color: blue
        color: white
        border-color: blue
        cursor: default
  .identifiers
    display: grid
    grid-template-columns: auto auto 1fr
    column-gap: 1rem
  .input
    .details
      display: grid
      grid-template-columns: auto 1fr
      column-gap: 1rem
      .label
        &:first-letter
          text-transform: uppercase

</style>