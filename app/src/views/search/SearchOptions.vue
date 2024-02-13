<!--
 Copyright 2023 Open Reaction Database Project Authors

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
-->

<script>
import ModalKetcher from '@/components/ModalKetcher'
import SearchItemList from './SearchItemList'
import MultiRangeSlider from "multi-range-slider-vue"

export default {
  components: {
    ModalKetcher,
    SearchItemList,
    MultiRangeSlider,
  },
  emits: ["searchOptions"],
  data() {
    return {
      test: null,
      showReagentOptions: false,
      showReactionOptions: false,
      showDatasetOptions: false,
      reagentOptions: {
        reactants: [],
        products: [],
        matchMode: "exact",
        useStereochemistry: false,
        similarityThreshold: 0.5,
      },
      reactionOptions: {
        reactionIds: [],
        reactionSmarts: [],
        min_yield: 50,
        max_yield: 100,
        min_conversion: 50,
        max_conversion: 100
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
      ketcherModalSet: "reactants",
      matchModes: ["exact", "similar", "substructure"]
    }
  },
  computed: {
    defaultQuery() {
      return this.$route.query
    },
    simThresholdDisplay() {
      let trailingZeros = ""
      const simThresh = this.reagentOptions.similarityThreshold.toString()
      if (simThresh.length < 2)
        trailingZeros = ".00"
      else if (simThresh.length < 4)
        trailingZeros = "0"
      return simThresh+trailingZeros
    }
  },
  methods: {
    openKetcherModal(componentSet, idx) {
      this.ketcherModalSmile = idx
      this.ketcherModalSet = componentSet
      this.showKetcherModal = true
    },
    emitSearchOptions() {
      const allComponents = [...this.reagentOptions.reactants, ...this.reagentOptions.products]
                              .map(component => ({...component, matchMode: this.reagentOptions.matchMode}))

      const searchOptions = {
        reagent: {
          reagents: allComponents,
          useStereochemistry: this.reagentOptions.useStereochemistry,
          similarityThreshold: this.reagentOptions.similarityThreshold
        },
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
        if (Array.isArray(q.component))
          q.component.forEach(comp => { this.addCompToOptions(comp) })
        else
          this.addCompToOptions(q.component)
        this.reagentOptions.useStereochemistry = q.use_stereochemistry || false
        this.reagentOptions.similarityThreshold = Number(q.similarity) || 0.5
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
      this.reactionOptions.min_yield = Number(q.min_yield) || 0
      this.reactionOptions.max_yield = Number(q.max_yield) || 100
      this.reactionOptions.min_conversion = Number(q.min_conversion) || 0
      this.reactionOptions.max_conversion = Number(q.max_conversion) || 100
      if (this.reactionOptions.reactionIds.length || this.reactionOptions.reactionSmarts.length || this.reactionOptions.yield !== 50 || this.reactionOptions.conversion !== 50) 
        this.showReactionOptions = true

      // general search params
      this.searchParams.limit = q.limit || 100
    },
    addCompToOptions (comp) {
      const compArray = comp.split(";")
      const compType = compArray[1] == "input" ? "reactants" : "products"
      this.reagentOptions.matchMode = compArray[2]
      this.reagentOptions[compType].push({smileSmart: compArray[0].replaceAll("%3D","="), source: compArray[1], matchMode: compArray[2]})
    },
    updateYield (e) {
      this.reactionOptions.min_yield = e.minValue
      this.reactionOptions.max_yield = e.maxValue
    },
    updateConversion (e) {
      this.reactionOptions.min_conversion = e.minValue
      this.reactionOptions.max_conversion = e.maxValue
    }
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
    span Components
    i.material-icons expand_less
  transition(name="expand")
    #searchByReagent.options-container(
      v-if='showReagentOptions'
    )
      .section
      .subtitle 
        .tabs
          .tab.capitalize(
            v-for='mode in matchModes'
            @click='reagentOptions.matchMode = mode'
            :class='reagentOptions.matchMode === mode ? "selected" : ""'
          ) {{mode}}
      .section
        .subtitle General Options
        .general.options
          label(for='stereo') Use Stereochemistry
          input#stereo(
            type='checkbox'
            v-model='reagentOptions.useStereochemistry'
          )
          template(v-if='reagentOptions.matchMode === "similar"')
            label(for='similarity') Similarity Threshold
            .slider-input
              .value {{simThresholdDisplay}}
              input#similarity(
                type='range'
                min="0.1"
                max="1.0"
                step="0.01"
                v-model='reagentOptions.similarityThreshold'
              )
      .section
        .subtitle Reactants & Reagents
        .reagent.options
          template(v-for='(reactant, idx) in reagentOptions.reactants')
            .draw
              button(@click='openKetcherModal("reactants", idx)')
                i.material-icons draw
            .field.long 
              input(
                type='text'
                v-model='reactant.smileSmart'
              )
            .delete
              button(@click='reagentOptions.reactants.splice(idx,1)')
                i.material-icons delete
          .copy(v-if='!reagentOptions.reactants.length') No components
          #add-component
            button(
              type='button' 
              @click='this.reagentOptions.reactants.push({smileSmart: "", source: "input"})'
            )
              i.material-icons add
              |  Add Component
      .section
        .subtitle Products
        .reagent.options
          template(v-for='(product, idx) in reagentOptions.products')
            .draw
              button(@click='openKetcherModal("products", idx)')
                i.material-icons draw
            .field.long 
              input(
                type='text'
                v-model='product.smileSmart'
              )
            .delete
              button(@click='reagentOptions.products.splice(idx,1)')
                i.material-icons delete
          .copy(v-if='!reagentOptions.products.length') No components
          #add-component
            button(
              type='button' 
              @click='this.reagentOptions.products.push({smileSmart: "", source: "output"})'
            )
              i.material-icons add
              |  Add Component
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
      .slider-input.multi
        label(for='yield') Yield
        .value {{reactionOptions.min_yield}}% - {{reactionOptions.max_yield}}%
        MultiRangeSlider(
          baseClassName="multi-range-slider"
          :miin='0'
          :max='100'
          :step='1'
          :ruler='false'
          :label='false'
          :minValue='reactionOptions.min_yield'
          :maxValue='reactionOptions.max_yield'
          @input='updateYield'
        )
      .slider-input.multi
        label(for='yield') Conversion
        .value {{reactionOptions.min_conversion}}% - {{reactionOptions.max_conversion}}%
        MultiRangeSlider(
          baseClassName="multi-range-slider"
          :miin='0'
          :max='100'
          :step='1'
          :ruler='false'
          :label='false'
          :minValue='reactionOptions.min_conversion'
          :maxValue='reactionOptions.max_conversion'
          @input='updateConversion'
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
  :smiles='reagentOptions[ketcherModalSet][ketcherModalSmile].smileSmart'
  @updateSmiles='(newSmiles) => reagentOptions[ketcherModalSet][ketcherModalSmile].smileSmart = newSmiles'
  @closeModal='showKetcherModal = false'
)
</template>

<style lang="sass" scoped>
@import '@/styles/vars'
@import '@/styles/transition'
@import '@/styles/tabs'
.search-options
  max-height: 90vh
  .options-title
    font-size: 1.5rem
    font-weight: 700
    cursor: pointer
    padding: 0.5rem 1rem
    border-top-left-radius: 0.25rem
    border-top-right-radius: 0.25rem
    margin-top: 1rem
    color: white
    background-color: $linkblue
    display: flex
    align-content: center
    justify-content: space-between
    width: 100%
    box-sizing: border-box
    transition: 0.25s
    i
      font-size: 2rem
      transition: 0.25s
    &.closed
      border-radius: 0.25rem
      i
        transform: rotate(180deg)
    &:first-child
      margin-top: 0
    &#searchParameters
      cursor: default
  .options-container
    background-color: white
    border-bottom-left-radius: 0.25rem
    border-bottom-right-radius: 0.25rem
    padding: 1rem
    margin-bottom: 1rem
    display: grid
    row-gap: 1rem
    border: 1px solid $medgrey
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
      grid-template-columns: auto 1fr auto
      column-gap: 1rem
      row-gap: 0.5rem
      align-items: center
      .label
        font-size: 1.1rem
      #add-component, .copy
        grid-column: 1 / 3
      .field.long
        min-width: 250px
        input
          width: 95%
    .general.options
      display: grid
      grid-template-columns: auto 1fr
      column-gap: 1rem
      row-gap: 0.5rem
      input
        max-width: 15px
        text-align: center
      .slider-input
        display: flex
        column-gap: 0.5rem
        #similarity
          // width: 8rem
          max-width: 8rem
  #searchByReaction
    .slider-input
      display: grid
      grid-template-columns: 10rem 2.5rem 1fr
      column-gap: 0.5rem
      .value
        text-align: right
      &.multi
        grid-template-columns: 4rem 6rem 1fr
        align-items: center
</style>