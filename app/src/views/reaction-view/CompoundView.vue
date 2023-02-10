<script>
import { reaction_pb } from "ord-schema"

export default {
  props: {
    component: Array,
  },
  watch: {
    compound: function (newVal, oldVal) {
      this.getCompoundSVG(newVal)
    }
  },
  data() {
    return {
      compoundSVG: null,
      compoundAmount: null,
    }
  },
  methods: {
    getCompoundSVG (component) {
      const compoundStr = component
                        .identifiersList
                        .find((identifier) => identifier.type == 2) //type 2 is SMILES
                        .value
      // prep compound
      const compound = new reaction_pb.Compound()
      const identifier = compound.addIdentifiers()
      identifier.setValue(compoundStr)
      identifier.setType(reaction_pb.CompoundIdentifier.CompoundIdentifierType.SMILES)
      // send http request
      return new Promise(resolve => {
        const xhr = new XMLHttpRequest()
        xhr.open("POST", "/api/render/compound/svg")
        const binary = compound.serializeBinary()
        xhr.responseType = "json"
        xhr.onload = function () {
          resolve(xhr.response)
        }
        xhr.send(binary)
      }).then(val => {
        this.compoundSVG = val
      })
    },
    getCompoundAmount (component) {
      console.log('amount',component.amount)
      const amount = new reaction_pb.Amount()
      console.log('amountconst',amount)
      const moles = amount
      return new Promise(resolve => {
        const xhr = new XMLHttpRequest()
        xhr.open("POST", "/api/render/compound/amount")
        const binary = amount.serializeBinary()
        xhr.responseType = "json"
        xhr.onload = function () {
          resolve(xhr.response)
        }
        xhr.send(binary)
      }).then(val => {
        this.compoundAmount = val
      })
    }
  },
  async mounted() {
    // console.log('component',this.component)
    this.getCompoundSVG(this.component[1].componentsList[0])  // TODO Can there be more than one component in here?
    this.getCompoundAmount(this.component[1].componentsList[0])
  }
}
</script>

<template lang="pug">
.details
  .svg(
    v-html='compoundSVG'
  )
  .amount {{compoundAmount}}
</template>

<style lang="sass" scoped>

</style>