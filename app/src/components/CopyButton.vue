<script>
export default {
  name: 'CopyButton',
  props: {
    textToCopy: String,
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
  button(@click='copy') Copy
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
  #copy-notification
    position: fixed
    background-color: $darkgrey
    color: white
    padding: 0.5rem 1rem
    border-radius: 0.25rem
</style>