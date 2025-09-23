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
import HeaderNav from "./components/HeaderNav"
import MainFooter from "./components/MainFooter"

export default {
  components: {
    HeaderNav,
    MainFooter,
  },
  computed: {
    noHeaderFooter() {
      const noHeaderRoutes = ['ketcher']
      return noHeaderRoutes.includes(this.$route.name)
    },
  }
}
</script>

<template lang="pug">
#main-container(
  :class='noHeaderFooter ? "full-height" : ""'
)
  HeaderNav(v-if='!noHeaderFooter')
  router-view
  MainFooter(v-if='!noHeaderFooter')
</template>

<style lang="sass">
@import url('https://fonts.googleapis.com/css2?family=Roboto:wght@400;700&display=swap')
@import '@/styles/vars'
body
  font-family: 'Roboto', sans-serif
  margin: 0
#app, #main-container
  min-height: 100vh
#main-container
  display: grid
  grid-template-rows: auto 1fr auto
  background-color: $lightgrey
  min-width: 800px
  &.full-height
    grid-template-rows: 1fr

button
  background-color: $linkblue
  border: none
  color: white
  border-radius: 0.25rem
  cursor: pointer
  padding: 0.25rem 0.5rem

// multi-value-slider styles
.multi-range-slider
  border: none !important
  box-shadow: none !important
  .thumb::before
    height: 18px !important
    width: 18px !important
    margin: -5px !important
    box-shadow: none !important
    background-color: #0d6efd !important
  .bar-left, .bar-inner, .bar-right
    box-shadow: none !important
    border: grey 1px solid
  .bar-inner
    background-color: #0d6efd !important
</style>
