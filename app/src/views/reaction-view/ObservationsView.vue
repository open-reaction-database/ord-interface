<script>
import { reaction_pb } from "ord-schema"

export default {
  props: {
    observations: Array,
  },
}
</script>

<template lang="pug">
.notes-view
  .details
    .label Time
    .label Comment
    template(v-for='obs in observations')
      .value {{obs.time}}
      .value {{obs.comment}}
</template>

<style lang="sass" scoped>
.notes-view
  .details
    display: grid
    grid-template-columns: auto 1fr
    column-gap: 1rem
    row-gap: 0.5rem
    .label
      font-weight: bold
</style>