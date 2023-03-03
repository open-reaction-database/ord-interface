<script>
import { reaction_pb } from "ord-schema"

export default {
  props: {
    events: Array,
  },
  data () {
    return {
      eventIdx: 0
    }
  },
  methods: {
  }
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