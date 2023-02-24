<script>
import { reaction_pb } from "ord-schema"
import CompoundView from "./CompoundView"
import FloatingModal from "../../components/FloatingModal"
import { amountObj, amountStr } from "../../utils/amount"

export default {
  props: {
    outcome: Object,
  },
  components: {
    CompoundView,
    FloatingModal,
  },
  data () {
    return {
      productsIdx: 0,
      showRawMeasurement: false,
      rawMeasurement: {},
    }
  },
  computed: {
    reactionTime () {
      const timeUnits = reaction_pb.Time.TimeUnit
      const type = Object.keys(timeUnits).find(key => timeUnits[key] == this.outcome.reactionTime.units)
      return `${this.outcome.reactionTime.value} ${type.toLowerCase()}${this.outcome.reactionTime.units !== 0 ? "(s)" : ""}`
    }
  },
  methods: {
    getMeasurementType (type) {
      const measurementTypes = reaction_pb.ProductMeasurement.ProductMeasurementType
      return Object.keys(measurementTypes).find(key => measurementTypes[key] == type)
    },
    getMeasurementValue (measurement) {
      if (measurement.percentage) 
        return `${measurement.percentage.value}%`
      if (measurement.amount) {
        console.log('amount',measurement.amount)
        return amountStr(amountObj(measurement.amount))
      }
      return ""
    },
    setRawMeasurement (measurement) {
      // make a deep copy of measurement and generate raw object 
      const raw = JSON.parse(JSON.stringify(measurement))
      raw.type = this.getMeasurementType(raw.type)
      if (raw.amount) {
        Object.keys(raw.amount).forEach(key => {
          if (raw.amount[key]) {
            raw.amount[key].units = amountObj(raw.amount).unitType
          }
        })
      }
      this.rawMeasurement = raw
      this.showRawMeasurement = true
    }
  }
}
</script>

<template lang="pug">
.outcomes-view
  template(v-if='outcome.reactionTime || outcome.conversion')
    .title Details
    .details
      template(v-if='outcome.reactionTime')
        .label Reaction Time
        .value {{reactionTime}}
      template(v-if='outcome.conversion')
        .label Conversion
        .value {{outcome.conversion}}
  .title Products
  .sub-section
    .tabs
      .tab(
        v-for='(product, idx) in outcome.productsList'
        @click='productsIdx = idx'
        :class='productsIdx === idx ? "selected" : ""'
      ) Product {{idx + 1}}
    .compound
      CompoundView(
        :component='outcome.productsList[productsIdx]'
      )
    .sub-title Measurements
    .measurements
      .label Type
      .label Value
      .label Analysis
      .label Raw
      template(v-for='measurement in outcome.productsList[productsIdx].measurementsList')
        .value {{getMeasurementType(measurement.type)}}
        .value#this-disappears {{getMeasurementValue(measurement)}}
        .value {{measurement.analysisKey}}
        .value 
          .raw 
            .button#this-button(@click='setRawMeasurement(measurement)') &lt;>
    floating-modal(
      v-if='showRawMeasurement'
      title="Raw Data"
      @closeModal='showRawMeasurement=false'
    )
      .data
        pre {{rawMeasurement}}

</template>

<style lang="sass" scoped>
@import "../../styles/vars"
@import "../../styles/tabs"
.outcomes-view
  .title
    font-weight: 700
    font-size: 1.5rem
    margin-bottom: 0.25rem
    margin-top: 1rem
  .details
    display: grid
    grid-template-columns: auto 1fr
    column-gap: 1rem
    row-gap: 0.5rem
  .sub-section
    border: 1px solid $medgrey
    border-radius: 0.25rem
    padding: 1rem
    .compound
      width: fit-content
    .sub-title
      font-size: 1.25rem
      font-weight: 700
    .measurements
      display: grid
      grid-template-columns: repeat(3, auto) 1fr
      column-gap: 0.5rem
      row-gap: 1rem
      padding: 0.5rem
      width: fit-content
      align-items: center
      .label
        font-weight: 700
      .raw
        .button
          padding: 0.5rem
          background-color: #a0a0a0
          color: white
          border-radius: 0.25rem
          cursor: pointer
</style>