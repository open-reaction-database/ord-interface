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
import LoadingSpinner from '@/components/LoadingSpinner'
import conditionUtil from '@/utils/conditions'
import outcomesUtil from '@/utils/outcomes'
import reaction_pb from 'ord-schema'
import CopyButton from '@/components/CopyButton'

export default {
  props: {
    reaction: Object,
    isSelected: Boolean,
    isSelectable: Boolean
  },
  components: {
    LoadingSpinner,
    CopyButton
  },
  data() {
    return {
      formattedResults: [],
      selectedReactions: [],
      reactionTable: null
    }
  },
  methods: {
    getReactionTable() {
      fetch(`/api/render/${this.reaction.reaction_id}`)
        .then(response => response.json())
        .then(responseData => {
          this.reactionTable = responseData
        })
    },
    getYield(measurements) {
      const yieldObj = measurements.find(m => m.type == 3) // ord-schema type 3 == "YIELD"
      if (yieldObj?.percentage) {
        return `${yieldObj.percentage.value}%`
      } else {
        return "Not listed"
      }
    },
    getConversion(reaction) {
      if (!reaction.outcomesList[0].conversion) return "Not listed"
      // decode conversion
    },
    conditionsAndDuration(reaction) {
      const details = []
      // get temp
      const tempSetPoint = conditionUtil.tempSetPoint(reaction.conditions.temperature?.setpoint)
      if (tempSetPoint !== "None")
        details.push(`at ${tempSetPoint}`)

      // get Pressure
      const pressureSetPoint = conditionUtil.pressureSetPoint(reaction.conditions.pressure?.setpoint)
      if (pressureSetPoint !== "None")
        details.push(`under ${pressureSetPoint}`)

      // get duration
      const formattedTime = outcomesUtil.formattedTime(reaction.outcomesList[0].reactionTime)
      if (formattedTime)
        details.push(`for ${formattedTime}`)

      return details
    },
    productIdentifier(identifier) {
      const identifierTypes = reaction_pb.CompoundIdentifier.CompoundIdentifierType
      const identifierType = Object.keys(identifierTypes).find(key => identifierTypes[key] == identifier.type)
      return `${identifierType}: ${identifier.value}`
    },
    updateSelectedReactions(event) {
      if (event.target.checked) {
        this.selectedReactions.push(event.target.value)
      } else {
        let idx = this.selectedReactions.indexOf(event.target.value)
        if (idx !== -1) {
          this.selectedReactions.splice(idx, 1)
        }
      }
    },
  },
  async mounted() {
    this.getReactionTable()
  },
}
</script>

<template lang="pug">
.reaction-container
  .reaction-row(:class='isSelected ? "selected" : ""')
    .select(v-if='isSelectable')
      input(
        type="checkbox"
        :id='"select_"+reaction.reaction_id'
        :value='reaction.reaction_id'
        :checked='isSelected'
        @change='$emit("clickedSelect", $event)'
      )
      label(:for='"select_"+reaction.reaction_id') Select reaction
    .reaction-table(
      v-html='reactionTable'
      v-if='reactionTable'
    )
    LoadingSpinner(v-else)
    .info
      .col.full
        router-link(
          :to='{ name: "reaction-view", params: {reactionId: reaction.reaction_id}}'
        ) 
          button View Full Details
      .col
        .yield Yield: {{getYield(reaction.data.outcomesList[0].productsList[0].measurementsList)}}
        .conversion Conversion: {{getConversion(reaction.data)}}
        .conditions Conditions: {{conditionsAndDuration(reaction.data).join("; ") || "Not Listed"}}
        .smile(v-if='reaction.data.outcomesList[0].productsList[0].identifiersList.length')
          CopyButton(
            :textToCopy='reaction.data.outcomesList[0].productsList[0].identifiersList[0].value'
          )
          .value Product {{productIdentifier(reaction.data.outcomesList[0].productsList[0].identifiersList[0])}}
      .col
        .creator Uploaded by {{reaction.data.provenance.recordCreated.person.name}}, {{reaction.data.provenance.recordCreated.person.organization}}
        .date Uploaded on {{new Date(reaction.data.provenance.recordCreated.time.value).toLocaleDateString()}}
        .doi DOI: {{reaction.data.provenance.doi}}
        .publication 
          a(
            :href='reaction.data.provenance.publicationUrl'
            target="_blank"
          ) Publication URL
        .dataset Dataset: {{reaction.dataset_id}}
</template>

<style lang="sass" scoped>
@import '@/styles/vars.sass'
.reaction-container
  text-decoration: none
  .reaction-row
    background-color: white
    border-radius: 0.25rem
    padding: 0.75rem
    margin-bottom: 1rem
    transition: 0.25s
    border: 4px solid white
    &:hover
      box-shadow: 0 0 5px $darkgrey
    .select
      input, label
        cursor: pointer
    .reaction-table
      color: black
      overflow-x: wrap
    .info
      display: grid
      grid-template-columns: repeat(2, 50%)
      row-gap: 0.5rem
      column-gap: 1rem
      margin-top: 1rem
      .col.full
        grid-column: 1/3
        button
          font-size: 1.2rem
      .col
        *
          margin-top: 0.5rem
      .smile
        display: flex
        column-gap: 0.5rem
        align-items: center
        margin-top: 0.25rem
        .value
          width: 100%
          white-space: nowrap
          overflow: hidden
          text-overflow: ellipsis
    &.selected
      border-color: $linkblue
</style>