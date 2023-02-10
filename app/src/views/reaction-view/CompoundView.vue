<script>
import { reaction_pb } from "ord-schema"

export default {
  props: {
    component: Object,
  },
  watch: {
    component: {
      handler(val) {
        this.getCompoundSVG(val)
      },
      deep: true,
      immediate: true,
    },
  },
  data() {
    return {
      compoundSVG: null,
    }
  },
  computed: {
    amountObj () {
      const amount = this.component.amount
      let units = {}
      let unitType = ""
      // determine unit type
      if (amount.moles) {
        units = reaction_pb.Moles.MolesUnit
        unitType="moles"
      } else if (amount.volume) {
        units = reaction_pb.Volume.VolumeUnit
        unitType = "volume"
      } else if (amount.mass) {
        units = reaction_pb.Mass.MassUnit
        unitType = "mass"
      } else if (amount.unmeasured) {
        console.log('unmeasured unit type')
        return {unitAmount: "", unitType: "unmeasured"}
      }
      const unitVal = amount[unitType].units
      const amountVal = amount[unitType].value
      const precision = amount[unitType].precision
      console.log('recision',precision)
      return {
        unitAmount: amountVal, 
        unitType: Object.keys(units).find(key => units[key] == unitVal)
      }
    },
    compoundAmount () {
      return `${Math.round(this.amountObj?.unitAmount * 1000) / 1000} ${this.amountObj?.unitType.toLowerCase()}`
    },
    compoundRole () {
      const role = this.component.reactionRole
      const types = reaction_pb.ReactionRole.ReactionRoleType
      return Object.keys(types).find(key => types[key] == role).toLowerCase()
    },
    rawData () {
      return {}
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
  },
  async mounted() {
    // console.log("schema",reaction_pb)
    // this.getCompoundSVG(this.component)
  }
}
</script>

<template lang="pug">
.compound-view
  .svg(
    v-html='compoundSVG'
  )
  .amount {{compoundAmount}}
  .role {{compoundRole}}
  .raw {{rawData}}
</template>

<style lang="sass" scoped>
.compound-view
  display: grid
  grid-template-columns: 1fr auto auto auto
  column-gap: 1rem
  align-items: center
</style>