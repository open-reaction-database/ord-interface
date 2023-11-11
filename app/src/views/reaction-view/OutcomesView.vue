<!--
 Copyright 2023 Open Reaction Database Project Authors

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
-->

<script>
import reaction_pb from "ord-schema"
import CompoundView from "./CompoundView"
import FloatingModal from "../../components/FloatingModal"
import { amountObj, amountStr } from "@/utils/amount"
import outcomesUtil from "@/utils/outcomes"

export default {
  props: {
    outcome: Object,
  },
  components: {
    CompoundView,
    FloatingModal
  },
  data () {
    return {
      productsIdx: 0,
      showRawMeasurement: false,
      rawMeasurement: {},
      analysesIdx: 0,
      showRawAnalysis: false,
      rawAnalysis: {}
    }
  },
  computed: {
    reactionTime () {
      return outcomesUtil.formattedTime(this.outcome.reactionTime)
    }
  },
  methods: {
    getMeasurementType (type) {
      const measurementTypes = reaction_pb.ProductMeasurement.ProductMeasurementType
      return Object.keys(measurementTypes).find(key => measurementTypes[key] == type)
    },
    getMeasurementValue (measurement) {
      if (measurement.percentage) 
        return `${Math.round(measurement.percentage.value * 10)/10}%`
      if (measurement.amount) {
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
    },
    getAnalysisType (type) {
      const analysisTypes = reaction_pb.Analysis.AnalysisType
      return Object.keys(analysisTypes).find(key => analysisTypes[key] == type)
    },
    setRawAnalysis (analysis) {
      const raw = JSON.parse(JSON.stringify(analysis))
      raw.type = this.getAnalysisType(raw.type)
      this.rawAnalysis = raw
      this.showRawAnalysis = true
    },
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
        .value {{outcome.conversion.value}} {{isNaN(outcome.conversion.precision) ? "" : `± ${outcome.conversion.precision}`}}
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
        .value {{getMeasurementValue(measurement)}}
        .value {{measurement.analysisKey}}
        .value 
          .raw 
            .button(@click='setRawMeasurement(measurement)') &lt;>
    floating-modal(
      v-if='showRawMeasurement'
      title="Raw Data"
      @closeModal='showRawMeasurement=false'
    )
      .data
        pre {{rawMeasurement}}
  template(v-if='outcome.analysesMap?.length')
    .title Analyses
    .sub-section
      .tabs
        .tab(
          v-for='(analysis, idx) in outcome.analysesMap'
          @click='analysesIdx = idx'
          :class='analysesIdx === idx ? "selected" : ""'
        ) {{analysis[0]}}
      .details
        .label Type
        .value {{getAnalysisType(outcome.analysesMap[analysesIdx][1].type)}}
        .label Details
        .value {{outcome.analysesMap[analysesIdx][1].details}}
        .label Raw
        .value
          .raw
            .button(@click='setRawAnalysis(outcome.analysesMap[analysesIdx][1])') &lt;>
    floating-modal(
      v-if='showRawAnalysis'
      title="Raw Data"
      @closeModal='showRawAnalysis=false'
    )
      .data
        pre {{rawAnalysis}}


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
    *
      display: flex
      align-items: center
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
        background-color: $darkgrey
        color: white
        border-radius: 0.25rem
        cursor: pointer
        width: fit-content
</style>