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
import reaction_pb from "ord-schema"

export default {
  props: {
    notes: Object,
  },
  computed: {
    notesToDisplay () {
      // grab fields with value and that aren't false for display
      return Object.keys(this.notes)
              .filter(key => this.notes[key])
              .map(key => {
                return {
                  val: this.notes[key],
                  label: this.camelToSpaces(key)
                }
              })
    }
  },
  methods: {
    camelToSpaces (str) {
      // convert camel case strings to string with spaces
      // also remove leading "is " where needed
      return str.replace(/([a-z])([A-Z])/g, '$1 $2')
              .replace("is ", "")
              .toLowerCase()
    }
  }
}
</script>

<template lang="pug">
.notes-view
  .details
    template(v-for='note in notesToDisplay')
      .label {{note.label}}
      .value {{note.val}}
</template>

<style lang="sass" scoped>
.notes-view
  .details
    display: grid
    grid-template-columns: auto 1fr
    column-gap: 1rem
    row-gap: 0.5rem
    .label:first-letter
      text-transform: capitalize
</style>