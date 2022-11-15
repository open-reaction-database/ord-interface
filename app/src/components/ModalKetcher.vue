<script>
export default {
  props: {
    smiles: String,
  },
  data() {
    return {
      contWin: null, //content window of iframe
      mutatedSmiles: "",
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
          }
        }
      }, 1000)
    },
    drawSmiles() {
      // this.ketcher.editor.struct(null);  // Clear any previous molecule.
      // const ketcherModal = document.getElementById('ketcher_modal');
      if (this.mutatedSmiles) {
        const xhr = new XMLHttpRequest();
        xhr.open('POST', '{{ url_for(".get_molfile") }}');
        xhr.responseType = 'json';
        xhr.onload = function () {
          if (xhr.status === 200) {
            const molblock = xhr.response;
            this.contWin.ketcher.setMolecule(molblock);
            // if (ketcherModal.hasClass('show')) {
            //   // If the modal is already open, we can simply set the molecule.
            //   this.ketcher.setMolecule(molblock);
            // } else {
            //   // Otherwise, we need to set up a callback, so that the molecule is set
            //   // only when Ketcher is open (to prevent a graphical glitch).
            //   ketcherModal.addEventListener('shown.bs.modal', function () {
            //     // This callback should only be ever run once, so make sure to remove it.
            //     ketcherModal.removeEventListener('shown.bs.modal');
            //     this.ketcher.setMolecule(molblock);
            //   });
            // }
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
    console.log('assigned',this.mutatedSmiles)
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
</template>

<style lang="sass" scoped>
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
        justify-content: end
        column-gap: 1rem
</style>