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
    },
    stirType () {
      const stirTypes = reaction_pb.StirringConditions.StirringMethodType
      return Object.keys(stirTypes).find(key => stirTypes[key] == this.conditions.stirring.type)
    },
    illumType () {
      const illum = this.conditions.illumination
      const illumTypes = reaction_pb.IlluminationConditions.IlluminationType
      const illumType = Object.keys(illumTypes).find(key => illumTypes[key] == illum.type)
      return `${illumType}${illum.details ? `: ${illum.details}` : ""}`
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

  .stirring.details(v-if='display==="stirring"')
    .label Type
    .value {{stirType}}
    template(v-if='conditions.stirring.details')
      .label Details
      .value {{conditions.stirring.details}}
      // TODO Flesh out stirring rate
    .label Rate
    .value {{conditions.stirring.rate || "UNSPECIFIED"}}
    template(v-if='conditions.stirring.rate?.rpm')
      .label RPM
      .value {{conditions.stirring.rate.rpm}}

  .illumination.details(v-if='display==="illumination"')
    .label Type
    .value {{illumType}}
    // TODO flesh out wave length
    .label Peak Wavelength
    .value {{conditions.illumination.peakWaveLength || "None"}}
    template(v-if='conditions.illumination.color')
      .label Color
      .value {{conditions.illumination.color}}
    // TODO flesh out distance
    .label Distance to Vessel
    .value {{conditions.illumination.distanceToVessel || "None"}}

  // TODO flesh out electrochemistry
  .electro.details(v-if='display === "electrochemistry"')
    .label Type
    .value {{conditions.electrochemistry}}
  
  // TODO flesh out flow
  .electro.details(v-if='display === "flow"')
    .label Type
    .value {{conditions.flow}}

  // TODO flesh out other
  .other.details(v-if='display === "other"')
    template(v-if='conditions.reflux')
      .label Reflux
      .value {{conditions.reflux}}
    template(v-if='conditions.ph')
      .label pH
      .value {{conditions.ph}}
    template(v-if='conditions.conditions_are_dynamic')
      .label Conditions are dynamic
      .value {{conditions.conditions_are_dynamic}}
    template(v-if='conditions.details')
      .label Details
      .value {{conditions.details}}
</template>

<style lang="sass" scoped>
.conditions-view
  .details
    display: grid
    grid-template-columns: auto 1fr
    column-gap: 1rem
    row-gap: 0.5rem
</style>