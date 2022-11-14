<script>
export default {
  methods: {
    drawSmiles() {
      // const ketcher = document.getElementById('ketcher-iframe').contentWindow.ketcher;
      // ketcher.editor.struct(null);  // Clear any previous molecule.
      // const ketcherModal = $('#ketcher_modal');
      // const smiles = $(node).val();
      // if (smiles) {
      //     const xhr = new XMLHttpRequest();
      //     xhr.open('POST', '{{ url_for(".get_molfile") }}');
      //     xhr.responseType = 'json';
      //     xhr.onload = function () {
      //         if (xhr.status === 200) {
      //             const molblock = xhr.response;
      //             if (ketcherModal.hasClass('show')) {
      //                 // If the modal is already open, we can simply set the molecule.
      //                 ketcher.setMolecule(molblock);
      //             } else {
      //                 // Otherwise, we need to set up a callback, so that the molecule is set
      //                 // only when Ketcher is open (to prevent a graphical glitch).
      //                 ketcherModal.on('shown.bs.modal', function () {
      //                     // This callback should only be ever run once, so make sure to remove it.
      //                     ketcherModal.off('shown.bs.modal');
      //                     ketcher.setMolecule(molblock);
      //                 });
      //             }
      //         }
      //     };
      //     xhr.send(smiles);
      // }
      // ketcher.successCallback = async function () {
      //     const smiles = await ketcher.getSmiles();
      //     if (smiles) {
      //         $(node).val(smiles);
      //     }
      // };
      // ketcherModal.modal('show');
    },
  },
  mounted() {
    console.log('show modal')
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
        ) Cancel
        button(
          type='button'
          onclick="document.getElementById('ketcher-iframe').contentWindow.ketcher.successCallback();"
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