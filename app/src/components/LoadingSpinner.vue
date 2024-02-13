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
  name: 'LoadingSpinner',
  data () {
    return {
    }
  },
}
</script>

<template lang="pug">
.spinner-main
  transition(
    name='fade'
  )
    .loader-ring
      div
      div
      div
      div
</template>

<style lang="sass" scoped>
@import '@/styles/vars.sass'
@import '@/styles/transition.sass'
.spinner-main
  width: 100%
  height: 100%
  display: flex
  .loader-ring
    display: inline-block
    position: relative
    width: 80px
    height: 80px
    margin: auto
    div 
      box-sizing: border-box
      display: block
      position: absolute
      width: 64px
      height: 64px
      margin: 8px
      border: 8px solid $linkblue
      border-radius: 50%
      animation: lds-ring 1.2s cubic-bezier(0.5, 0, 0.5, 1) infinite
      border-color: $linkblue transparent transparent transparent
      &:nth-child(1) 
        animation-delay: -0.45s
      &:nth-child(2) 
        animation-delay: -0.3s
      &:nth-child(3) 
        animation-delay: -0.15s

@keyframes lds-ring 
  0% 
    transform: rotate(0deg)
  100% 
    transform: rotate(360deg)
</style>