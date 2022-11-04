<script>
export default {
  data() {
    return {
      showReagentOptions: true,
      showReactionOptions: true,
      showDatasetOptions: true,
    }
  }
}
</script>

<template lang="pug">
.search-options
  .collapser(
    @click='showReagentOptions = !showReagentOptions'
    type='button'
    aria-controls='searchByReagent'
  ) Reagents
  #searchByReagent.options-container(
    v-if='showReagentOptions'
  )
    #reagents_content
      b(style='grid-column: text') SMILES/SMARTS
      b Source
      b Match Mode
      button#add_component(type='button' onclick='addComponent();' style='grid-column: text;')
        i.fas.fa-plus
        |  Add Component
      label(for='stereo' style='grid-column: text;') Use Stereochemistry
      input#stereo(type='checkbox')
      label(for='similarity' style='grid-column: text;') Similarity Threshold
      input#similarity(type='text' value='0.5' style='width: 60px;')
  .collapser(
    @click='showReactionOptions = !showReactionOptions'
    type='button' 
    aria-controls='searchByReaction'
  ) Reactions
  #searchByReaction.option-container(
    v-if='showReactionOptions'
  )
    #reactions_content
      label(for='reaction_ids') Reaction IDs
      textarea#reaction_ids
      label(for='reaction_smarts') Reaction SMARTS
      textarea#reaction_smarts
  .collapser(
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
  #go.mt-3
    label(for='limit') Result Limit
    input#limit(type='text' value='100' min='0' style='width: 60px')
</template>

<style lang="sass" scoped>
.search-options
  .collapser
    font-size: 1.5rem
    font-weight: 700
    cursor: pointer
</style>