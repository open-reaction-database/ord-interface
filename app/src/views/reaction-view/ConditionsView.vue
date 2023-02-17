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
      if (!setpoint) return "None"
      const tempUnits = reaction_pb.Temperature.TemperatureUnit
      const tempUnit = Object.keys(tempUnits).find(key => tempUnits[key] == setpoint.units)
      return `${setpoint.value} (± ${setpoint.precision}) ${tempUnit.toLowerCase()}`
    },
    pressureType () {
      const controlTypes = reaction_pb.PressureConditions.PressureControl.PressureControlType
      return Object.keys(controlTypes).find(key => controlTypes[key] == this.conditions.pressure.control.type)
    },
    pressureSetPoint () {
      const setpoint = this.conditions.pressure.setpoint
      if (!setpoint) return "None"
      const pressureUnits = reaction_pb.Pressure.PressureUnit
      const pressureUnit = Object.keys(pressureUnits).find(key => pressureUnits[key] == setpoint.units)
      return `${setpoint.value} (± ${setpoint.precision}) ${pressureUnit.toLowerCase()}`
    },
    pressureAtmo () {
      const atmo = this.conditions.pressure.atmosphere
      const atmoTypes = reaction_pb.PressureConditions.Atmosphere.AtmosphereType
      const atmoType = Object.keys(atmoTypes).find(key => atmoTypes[key] == atmo.type)
      return `${atmoType}${atmo.details ? `, ${atmo.details}` : ""}`
    }
  },
}
</script>

<template lang="pug">
.conditions-view
  .temperature.details(v-if='display==="temperature"')
    .label Control Type
    .value {{tempType}}
    template(v-if='conditions.temperature.control.details')
      .label Details
      .value {{conditions.temperature.control.details}}
    .label Setpoint
    .value {{tempSetPoint}}
    // TODO Flesh out temp measurements
    template(v-if='conditions.temperature.measurementsList?.length')
      .label Measurements
      .value {{conditions.temperature.measurementsList}}
  .pressure.details(v-if='display==="pressure"')
    .label Control Type
    .value {{pressureType}}
    template(v-if='conditions.pressure.control.details')
      .label Details
      .value {{conditions.pressure.control.details}}
    .label Setpoint
    .value {{pressureSetPoint}}
    .label Atmosphere
    .value {{pressureAtmo}}
    // TODO Flesh out pressure measurements
    template(v-if='conditions.pressure.measurementsList?.length')
      .label Measurements
      .value {{conditions.pressure.measurementsList}}
  //- .stirring.details(v-if='display==="stirring"')


</template>

<style lang="sass" scoped>
.conditions-view
  .details
    display: grid
    grid-template-columns: auto 1fr
    column-gap: 1rem
    row-gap: 0.5rem
</style>