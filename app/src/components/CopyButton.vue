<script>
export default {
  name: 'CopyButton',
  props: {
    textToCopy: String,
  },
  data () {
    return {
      displayNotification: false
    }
  },
  methods: {
    copy () {
      navigator.clipboard.writeText(this.textToCopy).then(() => {
        this.displayNotification = true
        setTimeout(() => {
          this.displayNotification = false
        }, 2000)
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
    #copy-notification Copied to clipboard!
</template>

<style lang="sass" scoped>
@import '@/styles/vars.sass'
@import '@/styles/transition.sass'
.copy-button-main
  #copy-notification
    position: fixed
    background-color: $darkgrey
    color: white
</style>