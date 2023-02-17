<script>
import { reaction_pb } from "ord-schema"
import FloatingModal from "../../components/FloatingModal"

export default {
  props: {
    component: Object,
  },
  components: {
    "floating-modal": FloatingModal
  },
  watch: {
    component: {
      handler(val) {
        this.getCompoundSVG(val)
        console.log('val',val)
      },
      deep: true,
      immediate: true,
    },
  },
  data() {
    return {
      compoundSVG: null,
      showRawData: false,
    }
  },
  computed: {
    amountObj () {
      const amount = this.component.amount
      let units = {}
      let unitCategory = ""
      // determine unit type
      if (amount.moles) {
        units = reaction_pb.Moles.MolesUnit
        unitCategory="moles"
      } else if (amount.volume) {
        units = reaction_pb.Volume.VolumeUnit
        unitCategory = "volume"
      } else if (amount.mass) {
        units = reaction_pb.Mass.MassUnit
        unitCategory = "mass"
      } else if (amount.unmeasured) {
        console.log('unmeasured unit type')
        return {unitAmount: "", unitCategory: "unmeasured"}
      }
      const unitVal = amount[unitCategory].units
      const amountVal = amount[unitCategory].value
      const precision = amount[unitCategory].precision
      console.log('recision',precision)
      return {
        unitAmount: amountVal, 
        unitType: Object.keys(units).find(key => units[key] == unitVal),
        unitCategory: unitCategory,
      }
    },
    compoundAmount () {
      return `${Math.round(this.amountObj?.unitAmount * 1000) / 1000} ${this.amountObj?.unitType.toLowerCase()}`
    },
    compoundRole () {
      const role = this.component.reactionRole
      const types = reaction_pb.ReactionRole.ReactionRoleType
      return Object.keys(types).find(key => types[key] == role)
    },
    rawData () {
      const returnObj = {}
      // set identifiers
      const idTypes = reaction_pb.CompoundIdentifier.CompoundIdentifierType
      returnObj["identifiers"] = this.component?.identifiersList?.map((identifier) => {
        return {
          type: Object.keys(idTypes).find(key => idTypes[key] == identifier.type),
          value: identifier.value
        }
      })
      // set amount
      returnObj["amount"] = {
        [this.amountObj.unitCategory]: {
          value: this.amountObj.unitAmount,
          units: this.amountObj.unitType,
        }
      }
      // set preparations
      if (this.component?.preparationsList?.length) {
        const prepTypes = reaction_pb.CompoundPreparation.CompoundPreparationType
        returnObj["preparations"] = this.component.preparationsList.map(prep => {
          return {
            type: Object.keys(prepTypes).find(key => prepTypes[key] == prep.type),
            details: prep.details,
          }
        })
      }
      returnObj["reaction_role"] = this.compoundRole
      return returnObj
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
  mounted() {
    // console.log("schema",reaction_pb)
  }
}
</script>

<template lang="pug">
.compound-view
  .svg(
    v-html='compoundSVG'
  )
  .amount {{compoundAmount}}
  .role {{compoundRole.toLowerCase()}}
  .raw 
    .button(@click='showRawData=true') &lt;> 
  floating-modal(
    v-if='showRawData'
    title="Raw Data"
    @closeModal='showRawData=false'
  )
    .data
      pre(v-for='identifier in rawData.identifiers') identifiers: {{identifier}}
      pre amount: {{rawData.amount}}
      pre reaction_role: {{rawData.reaction_role}}
      pre(v-for='prep in rawData.preparations') preparations: {{prep}}
</template>

<style lang="sass" scoped>
.compound-view
  display: grid
  grid-template-columns: 1fr auto auto auto
  column-gap: 1rem
  align-items: center
  .raw
    .button
      padding: 0.5rem
      background-color: #a0a0a0
      color: white
      border-radius: 0.25rem
      cursor: pointer
  .data
    width: 100%
    height: 100%
    overflow-x: auto
</style>