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
import outcomesUtil from "../../utils/outcomes"

export default {
  props: {
    observations: Array,
  },
  methods: {
    getTime(element) {
      return outcomesUtil.formattedTime(element)
    }
  },
}
</script>

<template lang="pug">
.notes-view
  .details
    .label Time
    .label Comment
    template(v-for='obs in observations')
      .value {{getTime(obs.time)}}
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