<script>
import { reaction_pb } from "ord-schema"

export default {
  props: {
    component: Object,
  },
  watch: {
    compound: function (newVal) {
      this.getCompoundSVG(newVal)
    }
  },
  data() {
    return {
      compoundSVG: null,
    }
  },
  computed: {
    amountUnit () {
      const amount = this.component.amount
      const molesUnits = reaction_pb.Moles.MolesUnit
      return Object.keys(molesUnits).find(key => molesUnits[key] == amount.moles.value)
    },
    compoundAmount () {
      return `${this.component.amount.moles.units} ${this.amountUnit.toLowerCase()}`
    },
    compoundRole () {
      const role = this.component.reactionRole
      const types = reaction_pb.ReactionRole.ReactionRoleType
      return Object.keys(types).find(key => types[key] == role).toLowerCase()
    },
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
  },
  async mounted() {
    console.log('schema',reaction_pb)
    this.getCompoundSVG(this.component)
  }
}
</script>

<template lang="pug">
.details
  .svg(
    v-html='compoundSVG'
  )
  .amount {{compoundAmount}}
  .role {{compoundRole}}
</template>

<style lang="sass" scoped>

</style>