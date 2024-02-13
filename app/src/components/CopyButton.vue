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
  name: 'CopyButton',
  props: {
    textToCopy: String,
    icon: {
      type: String,
      default: "content_copy"
    },
    buttonText: {
      type: String,
      default: ""
    }
  },
  data () {
    return {
      displayNotification: false,
      notificationStyle: {
        top: 0,
        left: 0
      }
    }
  },
  methods: {
    copy (event) {
      // move notification block to mouse cursor
      const { clientX, clientY } = event
      this.notificationStyle = {
        top: `${clientY}px`,
        left: `${clientX}px`
      }

      navigator.clipboard.writeText(this.textToCopy).then(() => {
        this.displayNotification = true

        setTimeout(() => {
          this.displayNotification = false
        }, 1500)
      })
    }
  }
}
</script>

<template lang="pug">
.copy-button-main
  button(@click='copy')
    i.material-icons {{icon}}
    .copy(v-if='buttonText') {{buttonText}}
  transition(
    name='fade'
  )
    #copy-notification(
      v-if='displayNotification'
      :style='notificationStyle'
    ) Copied to clipboard!
</template>

<style lang="sass" scoped>
@import '@/styles/vars.sass'
@import '@/styles/transition.sass'
.copy-button-main
  button
    padding: 0.25rem 0.5rem 
    display: flex
    column-gap: 0.25rem
    i
      font-size: 1rem
  #copy-notification
    position: fixed
    background-color: $darkgrey
    color: white
    padding: 0.5rem 1rem
    border-radius: 0.25rem
</style>