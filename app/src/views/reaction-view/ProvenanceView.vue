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
    provenance: Object,
  },
}
</script>

<template lang="pug">
.provenance-view
  template(v-if='provenance.experimenter')
    .title Experimenter
    .details.experimenter
      .label Username
      .label Name
      .label ORCID
      .label Organization
      .label Email
      .value {{provenance.experimenter.username}}
      .value {{provenance.experimenter.name}}
      .value {{provenance.experimenter.orcid}}
      .value {{provenance.experimenter.organization}}
      .value {{provenance.experimenter.email}}
  .details
    template(v-if='provenance.city')
      .label City
      .value {{provenance.city}}
    template(v-if='provenance.experimentStart')
      .label City
      .value {{new Date(provenance.experimentStart).toLocaleDateString()}}
    template(v-if='provenance.doi')
      .label DOI
      .value {{provenance.doi}}
    template(v-if='provenance.patent')
      .label Patent
      .value {{provenance.patent}}
    template(v-if='provenance.publicationUrl')
      .label Publication URL
      .value 
        a(:href='provenance.publicationUrl' target="_blank") {{provenance.publicationUrl}}
</template>

<style lang="sass" scoped>
.provenance-view
  .title
    font-weight: 700
    font-size: 1.5rem
    margin-bottom: 0.25rem
  .details
    display: grid
    grid-template-columns: auto 1fr
    column-gap: 1rem
    row-gap: 0.5rem
    &.experimenter
      grid-template-columns: repeat(4, auto) 1fr
      .label
        font-weight: 700
</style>