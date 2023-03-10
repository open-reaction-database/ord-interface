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
        reactionIds: [],
        reactionSmarts: [],
      },
      datasetOptions: {
        datasetIds: [],
        DOIs: []
      },
      searchParams: {
        limit: 100
      },
      showKetcherModal: false,
      ketcherModalSmile: 0,
    }
  },
  computed: {
    defaultQuery() {
      return this.$route.query
    },
  },
  methods: {
    openKetcherModal(idx) {
      this.ketcherModalSmile = idx
      this.showKetcherModal = true
    },
    emitSearchOptions() {
      const searchOptions = {
        reagent: this.reagentOptions,
        reaction: this.reactionOptions,
        dataset: this.datasetOptions,
        general: this.searchParams
      }
      this.$emit('searchOptions', searchOptions)
    },
    setDefaultValues() {
      const q = this.defaultQuery

      // reagent options
      if (q.component?.length) {
        q.component.forEach(comp => {
          const compArray = comp.split(";")
          this.reagentOptions.reagents.push({smileSmart: compArray[0], source: compArray[1], matchMode: compArray[2]})
        })
        this.reagentOptions.useStereochemistry = q.use_stereochemistry || false
        this.reagentOptions.similarityThreshold = q.similarity || 0.5
        this.showReagentOptions = true
      }

      // dataset options
      this.datasetOptions.datasetIds = q.dataset_ids?.split(",") || []
      this.datasetOptions.DOIs = q.dois?.split(",") || []
      if (this.datasetOptions.datasetIds.length || this.datasetOptions.DOIs.length) 
        this.showDatasetOptions = true

      // reaction options
      this.reactionOptions.reactionIds = q.reaction_ids?.split(",") || []
      this.reactionOptions.reactionSmarts = q.reaction_smarts?.split(",") || []
      if (this.reactionOptions.reactionIds.length || this.reactionOptions.reactionSmarts.length) 
        this.showReactionOptions = true

      // general search params
      this.searchParams.limit = q.limit || 100
    },
  },
  mounted() {
    this.setDefaultValues()
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
        :itemList.sync='reactionOptions.reactionIds'
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
        :itemList.sync='datasetOptions.datasetIds'
      )
      SearchItemList(
        title='DOIs'
        :itemList.sync='datasetOptions.DOIs'
      )
  #searchParameters.options-title Search Parameters
  .options-container
    .section
      label(for='limit') Result Limit 
      input#limit(
        type='number'
        min='0' 
        v-model='searchParams.limit'  
      )
    .search-button
      button(
        @click='emitSearchOptions'
      )
        b Search
  ModalKetcher(
    v-if='showKetcherModal'
    :smiles='reagentOptions.reagents[ketcherModalSmile].smileSmart'
    @updateSmiles='(newSmiles) => reagentOptions.reagents[ketcherModalSmile].smileSmart = newSmiles'
    @closeModal='showKetcherModal = false'
  )
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
    border-radius: 0.25rem
    padding: 1rem
    margin-bottom: 1rem
    display: grid
    row-gap: 1rem
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
    .search-button
      grid-column: 1 / 2
      margin-top: 1rem
      button
        font-size: 1.2rem
        padding: 0.5rem 1rem
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