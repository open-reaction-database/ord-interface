<script>
import { reaction_pb } from "ord-schema"

export default {
  props: {
    workup: Array,
  },
}
</script>

<template lang="pug">
.workups-view
  .details {{workup}}
</template>

<style lang="sass" scoped>
.workups-view
  .details
    display: grid
    grid-template-columns: auto 1fr
    column-gap: 1rem
    row-gap: 0.5rem
    .label
      font-weight: bold
</style>