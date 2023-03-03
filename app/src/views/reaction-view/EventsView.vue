<script>
import { reaction_pb } from "ord-schema"

export default {
  props: {
    events: Array,
  },
}
</script>

<template lang="pug">
.events-view Hello
  
</template>

<style lang="sass" scoped>
.events-view
  .details
    display: grid
    grid-template-columns: auto 1fr
    column-gap: 1rem
    row-gap: 0.5rem
</style>