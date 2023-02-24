<script>
import { reaction_pb } from "ord-schema"

export default {
  props: {
    outcome: Object,
  },
}
</script>

<template lang="pug">
.notes-view
  .title Details
  .details
    .label Hello
    .value {{outcome}}
</template>

<style lang="sass" scoped>
.notes-view
  .title
    font-weight: 700
    font-size: 1.5rem
    margin-bottom: 0.5rem
  .details
    display: grid
    grid-template-columns: auto 1fr
    column-gap: 1rem
    row-gap: 0.5rem
</style>