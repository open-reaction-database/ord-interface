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
#ketcher_modal.modal.fade(tabindex='-1' role='dialog' aria-hidden='true')
  .modal-dialog.modal-lg.modal-dialog-centered(role='document')
    .modal-content(style='width: fit-content;')
      .modal-body
        iframe#ketcher-iframe(src='/ketcher')
      .modal-footer
        button.btn.btn-secondary(type='button' data-bs-dismiss='modal') Cancel
        button.btn.btn-primary(type='button' data-bs-dismiss='modal' onclick="document.getElementById('ketcher-iframe').contentWindow.ketcher.successCallback();")
          | Save
</template>

<style lang="sass" scoped>

</style>