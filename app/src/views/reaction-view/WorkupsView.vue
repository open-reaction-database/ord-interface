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
    workup: Object,
  },
  computed: {
    workupType () {
      const workupTypes = reaction_pb.ReactionWorkup.ReactionWorkupType
      return Object.keys(workupTypes).find(key => workupTypes[key] == this.workup?.type)
    },
  }
}
</script>

<template lang="pug">
.workups-view
  .details
    .label Type
    .value {{workupType}}
    template(v-if='workup.details')
      .label Details
      .value {{workup.details}}
    template(v-if='workup.duration')
      .label Duration
      .value {{workup.duration}}
    template(v-if='workup.amount')
      .label Aliquot amount
      .value {{workup.amount}}
    template(v-if='workup.keepPhase')
      .label Phase kept
      .value {{workup.keepPhase}}
    template(v-if='workup.targetPh')
      .label Target pH
      .value {{workup.targetPh}}
    template(v-if='workup.isAutomated')
      .label Automated
      .value {{workup.Automated}}
</template>

<style lang="sass" scoped>
.workups-view
  .details
    display: grid
    grid-template-columns: auto 1fr
    column-gap: 1rem
    row-gap: 0.5rem
</style>