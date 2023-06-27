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
  },
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
    .value {{vessel.material || "UNSPECIFIED"}}
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
    .value {{setup.environment?.type || "UNSPECIFIED"}}
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