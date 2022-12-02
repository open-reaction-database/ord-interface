<script>
export default {
  data() {
    return {
      reaction: {},
      loading: true
    }
  },
  computed: {
    reactionId() {
      return this.$route.params.reactionId
    }
  },
  methods: {
    getReactionData () {
      fetch(`/api/id/${this.reactionId}`, {method: "GET"})
        .then(response => response.json())
        .then(data => {
          this.reaction = data
          this.loading = false
          console.log('reaction',this.reaction)
        })
    }
  },
  mounted() {
    this.getReactionData()
  }
}
</script>

<template lang="pug">
.main-reaction-view Hello World
</template>

<style lang="sass" scoped>

</style>