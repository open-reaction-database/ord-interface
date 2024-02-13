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
import conditionUtil from "@/utils/conditions"

export default {
  props: {
    conditions: Object,
    display: String,
  },
  computed: {
    tempType () {
      return conditionUtil.tempType(this.conditions.temperature.control.type)
    },
    tempSetPoint () {
      return conditionUtil.tempSetPoint(this.conditions.temperature.setpoint)
    },
    pressureType () {
      return conditionUtil.pressureType(this.conditions.pressure.control.type)
    },
    pressureSetPoint () {
      return conditionUtil.pressureSetPoint(this.conditions.pressure?.setpoint)
    },
    pressureAtmo () {
      return conditionUtil.pressureAtmo(this.conditions.pressure.atmosphere)
    },
    stirType () {
      return conditionUtil.stirType(this.conditions.stirring.type)
    },
    stirRate () {
      return conditionUtil.stirRate(this.conditions.stirring.rate)
    },
    illumType () {
      return conditionUtil.illumType(this.conditions.illumination)
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
    .label Rate
    .value {{stirRate || "UNSPECIFIED"}}
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