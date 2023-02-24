<script>
import { reaction_pb } from "ord-schema"

export default {
  props: {
    outcome: Object,
  },
  data () {
    return {
      productsIdx: 0,
    }
  },
  computed: {
    reactionTime () {
      const timeUnits = reaction_pb.Time.TimeUnit
      const type = Object.keys(timeUnits).find(key => timeUnits[key] == this.outcome.reactionTime.units)
      return `${this.outcome.reactionTime.value} ${type.toLowerCase()}${this.outcome.reactionTime.units !== 0 ? "(s)" : ""}`
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
</style>