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
  width: 100%
  height: 100%
  .modal-background
    position: absolute
    width: 100%
    height: 100%
    top: 0
    left: 0
    background-color: #00000088
  .modal-holder
    position: relative
    height: 50%
    pointer-events: none
    .modal-container
      position: absolute
      max-width: 80%
      min-width: 50%
      height: 60vh
      background-color: white
      border-radius: 0.25rem
      margin: auto
      padding: 1rem 0
      top: 40%
      -ms-transform: translateY(-40%)
      transform: translateY(-40%)
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
        overflow-y: auto
        height: calc(100% - 3.5rem)
</style>