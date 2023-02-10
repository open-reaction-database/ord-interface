<script>
import { reaction_pb } from "ord-schema"

export default {
  props: {
    compound: String,
  },
  data() {
    return {
      compoundSVG: "test"
    }
  },
  computed: {
  },
  methods: {
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
          resolve(xhr.response)
        }
        xhr.send(binary)
      }).then(val => {
        this.compoundSVG = val
      })
      // return compound
    }
  },
  async mounted() {
    console.log('compound',this.compound)
    this.getCompoundSVG(this.compound)
  }
}
</script>

<template lang="pug">
.details
  .svg(
    v-html='compoundSVG'
  )
</template>

<style lang="sass" scoped>

</style>