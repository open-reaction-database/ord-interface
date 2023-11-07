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
    events: Array,
  },
  data () {
    return {
      eventIdx: 0
    }
  },
}
</script>

<template lang="pug">
.events-view
  .tabs
    .tab(
      v-for='(event, idx) in events'
      @click='eventIdx = idx'
      :class='eventIdx === idx ? "selected" : ""'
    ) {{new Date(event.time.value).toLocaleString()}}
  .details
    template(v-if='events[eventIdx].details')
      .label Details
      .value {{events[eventIdx].details}}
    template(v-if='events[eventIdx].person?.username')
      .label Username
      .value {{events[eventIdx].person.username}}
    template(v-if='events[eventIdx].person?.name')
      .label Name
      .value {{events[eventIdx].person.name}}
    template(v-if='events[eventIdx].person?.orcid')
      .label ORCID
      .value {{events[eventIdx].person.orcid}}
    template(v-if='events[eventIdx].person?.organization')
      .label Organization
      .value {{events[eventIdx].person.organization}}
    template(v-if='events[eventIdx].person?.email')
      .label Email
      .value {{events[eventIdx].person.email}}
  
</template>

<style lang="sass" scoped>
@import "../../styles/tabs"
.events-view
  .details
    display: grid
    grid-template-columns: auto 1fr
    column-gap: 1rem
    row-gap: 0.5rem
</style>