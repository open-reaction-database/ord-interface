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

export default {
  props: {
    setup: Object,
    display: String,
  },
  computed: {
    vessel () {
      return this.setup.vessel
    },
    vesselType () {
      const vesselTypes = reaction_pb.Vessel.VesselType
      return Object.keys(vesselTypes).find(key => vesselTypes[key] == this.vessel.type)
    },
    vesselVolume () {
      const unitLabels = reaction_pb.Volume.VolumeUnit
      const label = Object.keys(unitLabels).find(key => unitLabels[key] == this.vessel.volume.units)
      return `${this.vessel.volume.value} ${label.toLowerCase()}`
    },
    vesselAttachments () {
      const attachTypes = reaction_pb.VesselAttachment.VesselAttachmentType
      return this.vessel.attachmentsList.map(attach => {
        const type = Object.keys(attachTypes).find(key => attachTypes[key] == attach.type)
        return `${type}${attach.details ? `: ${attach.details}` : ""}`
      }).join(", ")
    },
    vesselPrep () {
      const prepTypes = reaction_pb.VesselPreparation.VesselPreparationType
      return this.vessel.preparationsList.map(prep => {
        const type = Object.keys(prepTypes).find(key => prepTypes[key] == prep.type)
        return `${type}${prep.details ? `: ${prep.details}` : ""}`
      }).join(", ")
    },
    vesselMaterial () {
      const materialTypes = reaction_pb.VesselMaterial.VesselMaterialType
      return Object.keys(materialTypes).find(key => materialTypes[key] == this.vessel.material?.type)
    },
    reactionEnv () {
      const envVal = this.setup.environment?.type
      if (!envVal) return null
      const envTypes = reaction_pb.ReactionSetup.ReactionEnvironment.ReactionEnvironmentType
      return Object.keys(envTypes).find(key => envTypes[key] == envVal)
    }
  },
  mounted () {
    console.log('pb',reaction_pb)
    console.log('props',this.setup)
  }
}
</script>

<template lang="pug">
.setup-view
  .vessel.details(v-if='display=="vessel"')
    .label Type
    .value {{vesselType}}
    template(v-if='vessel.details')
      .label Details
      .value {{vessel.details}}
    .label Material
    .value {{vesselMaterial || "UNSPECIFIED"}}
    .label Volume
    .value {{vesselVolume}}
    template(v-if='vessel.attachmentsList?.length')
      .label Attachments
      .value {{vesselAttachments}}
    template(v-if='vessel.preparationsList?.length')
      .label Preparations
      .value {{vesselPrep}}
  .environment.details(v-if='display=="environment"')
    .label Type
    .value {{reactionEnv || "UNSPECIFIED"}}
    template(v-if='setup.environment?.details')
      .label Details
      .value {{setup.environment.details}}
  .automated.details(v-if='display=="automated"')
    .label Platform
    .value {{setup.automation_platform || "UNSPECIFIED"}}
</template>

<style lang="sass" scoped>
.setup-view
  .details
    display: grid
    grid-template-columns: auto 1fr
    column-gap: 1rem
    row-gap: 0.5rem
</style>