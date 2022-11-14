<script>
import ModalKetcher from '@/components/ModalKetcher'
import SearchItemList from './SearchItemList'

export default {
  components: {
    ModalKetcher,
    SearchItemList
  },
  data() {
    return {
      showReagentOptions: false,
      showReactionOptions: false,
      showDatasetOptions: false,
      reagentOptions: {
        reagents: [],
        useStereochemistry: false,
        similarityThreshold: 0.5,
      },
      reactionOptions: {
        reactionIDs: [],
        reactionSmarts: [],
      },
      datasetOptions: {
        datasetIDs: [],
        dois: []
      },
      searchParams: {
        limit: 100
      },
      showKetcherModal: true,
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
    :class='showReagentOptions ? "" : "closed"'
  ) 
    span Reagents
    i.material-icons expand_less
  transition(name="expand")
    #searchByReagent.options-container(
      v-if='showReagentOptions'
    )
      .section
        .subtitle Components
        .reagent.options
          .spacer
          .label SMILES/SMARTS
          .label Source
          .label Match Mode
          .spacer
          template(v-for='(reagent, idx) in reagentOptions.reagents')
            .draw
              button(@click='openKetcherModal(idx)')
                i.material-icons draw
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
          .copy(v-if='!reagentOptions.reagents.length') No components
          #add-component
            button(
              type='button' 
              @click='this.reagentOptions.reagents.push({smileSmart: "", source: "input", matchMode: "exact"})'
            )
              i.material-icons add
              |  Add Component
      .section
        .subtitle General Options
        .general.options
          label(for='stereo') Use Stereochemistry
          input#stereo(s
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
    :class='showReactionOptions ? "" : "closed"'
  ) 
    span Reactions
    i.material-icons expand_less
  transition(name="expand")
    #searchByReaction.options-container(
      v-if='showReactionOptions'
    )
      SearchItemList(
        title='Reaction IDs'
        :itemList.sync='reactionOptions.reactionIDs'
      )
      SearchItemList(
        title='Reaction SMARTS'
        :itemList.sync='reactionOptions.reactionSmarts'
      )
  .options-title(
    @click='showDatasetOptions = !showDatasetOptions'
    :class='showDatasetOptions ? "" : "closed"'
  ) 
    span Datasets
    i.material-icons expand_less
  transition(name="expand")
    #searchByDataset.options-container(
      v-if='showDatasetOptions'
    )
      SearchItemList(
        title='Dataset IDs'
        :itemList.sync='reactionOptions.datasetIDs'
      )
      SearchItemList(
        title='DOIs'
        :itemList.sync='reactionOptions.DOIs'
      )
  #searchParamaters.options-title Search Paramaters
  .options-container
    .section
      label(for='limit') Result Limit 
      input#limit(
        type='number'
        min='0' 
        v-model='searchParams.limit'  
      )
  ModalKetcher(v-if='showKetcherModal')
</template>

<style lang="sass" scoped>
@import '@/styles/vars'
@import '@/styles/transition'
.search-options
  .options-title
    font-size: 1.5rem
    font-weight: 700
    cursor: pointer
    padding: 0.5rem 1rem
    border-radius: 0.25rem
    margin-bottom: 1rem
    color: white
    background-color: $linkblue
    display: flex
    align-content: center
    justify-content: space-between
    width: 13rem
    i
      font-size: 2rem
      transition: 0.25s
    &.closed
      i
        transform: rotate(180deg)
  .options-container
    background-color: white
    width: 100%
    border-radius: 0.25rem
    padding: 1rem
    margin-bottom: 1rem
    display: grid
    grid-template-columns: 1fr 1fr
    .subtitle
      font-size: 1.25rem
      font-weight: 700
      margin-bottom: 0.5rem
    .options
      margin-left: 1rem
    button
      display: flex
      align-items: center
      i
        font-size: 1.1rem
    input, select
      font-size: 1rem
  #searchByReagent
    .reagent.options
      display: grid
      grid-template-columns: repeat(4, auto) 1fr
      column-gap: 1rem
      row-gap: 0.5rem
      align-items: center
      .label
        font-size: 1.1rem
      #add-component, .copy
        grid-column: 1 / 3
    .general.options
      display: grid
      grid-template-columns: auto 1fr
      column-gap: 1rem
      row-gap: 0.5rem
      input
        max-width: 50px
        text-align: center
</style>