<script>
import KetcherModal from '@/components/KetcherModal'

export default {
  components: {
    KetcherModal
  },
  data() {
    return {
      showReagentOptions: true,
      showReactionOptions: true,
      showDatasetOptions: true,
      reagentOptions: {
        reagents: [{smileSmart: "", source: "input", matchMode: "exact"}],
        useStereochemistry: false,
        similarityThreshold: 0.5,
      },
      showKetcherModal: false,
    }
  },
  methods: {
    openKetcherModal(idx) {
      this.showKetcherModal = true
    }
  },
  mounted() {
  }
}
</script>

<template lang="pug">
.search-options
  .options-title(
    @click='showReagentOptions = !showReagentOptions'
    type='button'
    aria-controls='searchByReagent'
  ) Reagents
  #searchByReagent.options-container(
    v-if='showReagentOptions'
  )
    .reagent-options
      .spacer
      .label SMILES/SMARTS
      .label Source
      .label Match Mode
      .spacer
      template(v-for='(reagent, idx) in reagentOptions.reagents')
        .draw
          button(@click='openKetcherModal(idx)')
            i.material-icons edit
        .field.long 
          input(
            type='text'
            v-model='reagent.smileSmart'
          )
        .field
          select(v-model='reagent.source')
            option(value='input') input
            option(value='output') output
        .field
          select(v-model='reagent.matchMode')
            option(value='exact') exact
            option(value='similar') similar
            option(value='substructure') substructure
            option(value='smarts') smarts
        .delete
          button(@click='reagentOptions.reagents.splice(idx,1)')
            i.material-icons delete
    button#add_component(
      type='button' 
      @click='this.reagentOptions.reagents.push({smileSmart: "", source: "input", matchMode: "exact"})'
    )
      i.material-icons add
      |  Add Component
    .general-options
      label(for='stereo') Use Stereochemistry
      input#stereo(
        type='checkbox'
        v-model='reagentOptions.useStereochemistry'
      )
      label(for='similarity') Similarity Threshold
      input#similarity(
        type='number'
        step=0.1
        v-model='reagentOptions.similarityThreshold'
        )
  .options-title(
    @click='showReactionOptions = !showReactionOptions'
    type='button' 
    aria-controls='searchByReaction'
  ) Reactions
  #searchByReaction.options-container(
    v-if='showReactionOptions'
  )
    #reactions_content
      label(for='reaction_ids') Reaction IDs
      textarea#reaction_ids
      label(for='reaction_smarts') Reaction SMARTS
      textarea#reaction_smarts
  .options-title(
    @click='showDatasetOptions = !showDatasetOptions'
    type='button' 
    aria-controls='searchByDataset'
  ) Datasets
  #searchByDataset.options-container(
    v-if='showDatasetOptions'
  )
    #datasets_content
      label(for='dataset_ids') Dataset IDs
      textarea#dataset_ids
      label(for='dois') DOIs
      textarea#dois
  #searchParamaters.options-title Search Paramaters
  .options-container
    label(for='limit') Result Limit
    input#limit(type='text' value='100' min='0' style='width: 60px')
  KetcherModal(v-if='showKetcherModal')
</template>

<style lang="sass" scoped>
.search-options
  .options-title
    font-size: 1.5rem
    font-weight: 700
    cursor: pointer
    margin-bottom: 1rem
  .options-container
    background-color: white
    width: 100%
    border-radius: 0.25rem
    padding: 1rem
    margin-bottom: 1rem
    button
      display: flex
      align-items: center
      i
        font-size: 1.1rem
    input, select
      font-size: 1rem
  #searchByReagent
    .reagent-options
      display: grid
      grid-template-columns: repeat(4, auto) 1fr
      column-gap: 1rem
      row-gap: 0.5rem
      margin-bottom: 0.5rem
      align-items: center
      .label
        font-weight: 700
    .general-options
      display: grid
      grid-template-columns: auto 1fr
      column-gap: 1rem
      row-gap: 0.5rem
      margin-top: 1rem
      input
        max-width: 50px
        text-align: center
</style>