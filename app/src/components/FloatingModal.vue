<script>
export default {
  name: "FloatingModal",
  props: {
    title: String,
  },
  data() {
    return {};
  },
  methods: {
    closeModal() {
      this.$emit("closeModal");
    },
  },
};
</script>

<template lang="pug">
.modal-main
  .modal-background(@click='closeModal()')
  .modal-holder
    .modal-container
      .header 
        .title {{title}}
        .close(@click='closeModal()') &#10005;
      .body
        slot
</template>

<style lang="sass" scoped>
@import '@/styles/transition.sass'
.modal-main
  position: fixed
  top: 0
  left: 0
  width: 100vw
  height: 100vh
  .modal-background
    position: absolute
    width: 100vw
    height: 100vh
    top: 0
    left: 0
    background-color: #00000088
  .modal-holder
    position: relative
    height: 50%
    pointer-events: none
    .modal-container
      position: absolute
      max-width: 80vw
      min-width: 50vw
      max-height: 80vh
      min-height: 50vh
      background-color: white
      border-radius: 0.25rem
      margin: auto
      padding: 1rem 0
      top: 50%
      -ms-transform: translateY(-50%)
      transform: translateY(-50%)
      left: 50%
      -ms-transform: translateX(-50%)
      transform: translateX(-50%)
      pointer-events: auto
      .header
        font-weight: 700
        font-size: 1.5rem
        display: grid
        grid-template-columns: 1fr auto
        border-bottom: 2px solid black
        padding: 0 1rem 1rem
        .close
          cursor: pointer
      .body
        padding: 1rem 1rem 0
</style>