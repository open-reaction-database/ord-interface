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
    getReactionData (reactionId) {
      // fetch(`/api/id/${this.reactionId}`, {method: "GET"})
      //   .then(response => response.json())
      //   .then(data => {
      //     this.reaction = data
      //     this.loading = false
      //     console.log('reaction',this.reaction)
      //   })
      return new Promise(resolve => {
        const xhr = new XMLHttpRequest();
        xhr.open('GET', '/editor/reaction/id/' + reactionId + '/proto');
        xhr.responseType = 'arraybuffer';
        xhr.onload = function() {
          console.log(xhr.response)
          // asserts.assertInstanceof(xhr.response, ArrayBuffer); // Type hint.
          // const bytes = new Uint8Array(xhr.response);
          // const reaction = Reaction.deserializeBinary(bytes);
          // resolve(reaction);
        };
        xhr.send();
      });
    },
  },
  async mounted() {
    this.reaction = await this.getReactionData(this.reactionId)
  }
}
</script>

<template lang="pug">
.main-reaction-view Hello World
</template>

<style lang="sass" scoped>

</style>