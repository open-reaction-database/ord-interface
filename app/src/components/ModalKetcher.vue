<script>
import LoadingSpinner from '@/components/LoadingSpinner'

export default {
  components: {
    'LoadingSpinner': LoadingSpinner,
  },
  props: {
    smiles: String,
  },
  data() {
    return {
      contWin: null, //content window of iframe
      mutatedSmiles: "",
      loading: true,
    }
  },
  methods: {
    getKetcher() {
      const getKetcherInterval = setInterval(() => {
        //attempt to get Ketcher from the iframe once it's loaded
        if (!this.contWin) {
          if (document.getElementById('ketcher-iframe').contentWindow.ketcher) {
            // assigning contentWindow.ketcher doesn't seem to persist as expected
            // so we just assign the contentWindow and reference .ketcher
            this.contWin = document.getElementById('ketcher-iframe').contentWindow
            clearInterval(getKetcherInterval)
            this.drawSmiles()
            this.loading = false
          }
        }
      }, 1000)
    },
    drawSmiles() {
      if (this.mutatedSmiles) {
        // get molblock if we already have SMILES
        const xhr = new XMLHttpRequest();
        xhr.open('POST', 'api/ketcher/molfile');
        xhr.responseType = 'json';
        xhr.onload = function () {
          if (xhr.status === 200) {
            const molblock = xhr.response;
            document.getElementById('ketcher-iframe').contentWindow.ketcher.setMolecule(molblock);
          }
        };
        xhr.send(this.mutatedSmiles);
      }
    },
    async saveSmiles() {
      this.mutatedSmiles = await this.contWin.ketcher.getSmiles();
      this.$emit('updateSmiles', this.mutatedSmiles)
      this.$emit('closeModal')
    },
  },
  mounted() {
    this.mutatedSmiles = this.smiles
    this.getKetcher()
  }
}
</script>

<template lang="pug">
.background
  #ketcher_modal.modal
    .modal-content
      .modal-body
        iframe#ketcher-iframe( src='/ketcher' )
      .modal-footer
        button(
          type='button'
          @click='this.$emit("closeModal")'
        ) Cancel
        button(
          type='button'
          @click='saveSmiles()'
        ) Save
      transition(name="fade")
        .modal-loading(
          v-if='loading'
        )
          LoadingSpinner
</template>

<style lang="sass" scoped>
@import '@/styles/transition.sass'
.background
  width: 100vw
  height: 100vh
  z-index: 1
  background-color: #00000099
  position: absolute
  top: 0
  left: 0
  .modal
    width: 75vw
    height: 75vh
    background-color: white
    border-radius: 0.25rem
    position: absolute
    top: 50%
    left: 50%
    transform: translate(-50%, -50%)
    overflow: hidden
    max-width: 1000px
    max-height: 800px
    .modal-content
      width: 100%
      height: 100%
      display: grid
      grid-template-rows: 1fr auto
      .modal-body
        width: 100%
        height: 100%
        iframe
          width: 100%
          height: 100%
          border: 0
      .modal-footer
        padding: 1rem
        border-top: 2px solid #eff2f5 // matches ketcher bottom bar
        display: flex
        justify-content: flex-end
        column-gap: 1rem
    .modal-loading
      background-color: white
      width: 100%
      height: 100%
      position: absolute
</style>