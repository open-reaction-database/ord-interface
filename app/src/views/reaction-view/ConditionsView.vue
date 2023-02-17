<script>
import { reaction_pb } from "ord-schema"

export default {
  props: {
    conditions: Object,
    display: String,
  },
  computed: {
    tempType () {
      const controlTypes = reaction_pb.TemperatureConditions.TemperatureControl.TemperatureControlType
      return Object.keys(controlTypes).find(key => controlTypes[key] == this.conditions.temperature.control.type)
    },
    tempSetPoint () {
      const setpoint = this.conditions.temperature.setpoint
      return `${setpoint.value} (± ${setpoint.precision}) ${setpoint.units}`
    }
  },
}
</script>

<template lang="pug">
.conditions-view
  .temperature.details(v-if='display=="temperature"')
    .label Control Type
    .value {{tempType}}
    .label Setpoint
    .value {{tempSetPoint}}
</template>

<style lang="sass" scoped>
.conditions-view
  .details
    display: grid
    grid-template-columns: auto 1fr
    column-gap: 1rem
    row-gap: 0.5rem
</style>