<script>
import { reaction_pb } from "ord-schema"
import CompoundView from "./CompoundView"
import SetupView from "./SetupView"
import ConditionsView from "./ConditionsView"
import NotesView from "./NotesView"
import ObservationsView from "./ObservationsView"
import WorkupsView from "./WorkupsView"

export default {
  components: {
    CompoundView,
    SetupView,
    ConditionsView,
    NotesView,
    ObservationsView,
    WorkupsView,
  },
  data() {
    return {
      reaction: {},
      reactionSummary: null,
      reactionBytes: null,
      loading: true,
      inputsIdx: 0,
      setupTabs: [
        "vessel",
        "environment",
        "automation",
      ],
      setupTab: "vessel",
      conditionTabs: [
        "temperature",
        "pressure",
        "stirring",
        "illumination",
        "electrochemistry",
        "flow",
        "other",
      ],
      conditionTab: "temperature",
      workupsTab: 0,
    }
  },
  computed: {
    reactionId() {
      return this.$route.params.reactionId
    },
    displayInputs() {
      if (!this.reaction) return {}
      let returnArr = {...this.reaction.inputsMap[this.inputsIdx][1]}
      // filter out null/undefined values and arrays
      return Object.fromEntries(Object.entries(returnArr).filter(([_,v]) => v != null && !Array.isArray(v)))
    },
    displayConditionsOther() {
      const otherFields = [
        "reflux",
        "ph",
        "conditions_are_dynamic",
        "details",
      ]
      return otherFields.find(key => this.reaction.conditions[key])
    }
  },
  methods: {
    getReactionData () {
      return new Promise(resolve => {
        const xhr = new XMLHttpRequest();
        xhr.open("GET", `/api/getReaction/${this.reactionId}`)
        xhr.responseType = "arraybuffer";
        xhr.onload = () => {
          // if response is good, deserialize reaction and return object from protobuff
          let reaction = null
          if (xhr.response !== null) {
            this.reactionBytes = new Uint8Array(xhr.response);
            reaction = reaction_pb.Reaction.deserializeBinary(this.reactionBytes).toObject();
            // sort inputs by addition order
            reaction.inputsMap.sort((a,b) => a[1].additionOrder - b[1].additionOrder)
          }
          resolve(reaction);
        }
        xhr.send()
      })
    },
    async getReactionSummary () {
      const res = await fetch(`/api/render/${this.reactionId}?compact=false`)
      const data = await res.json()
      return data
    },
    getReactionType (id) {
      const identifiers = reaction_pb.ReactionIdentifier.ReactionIdentifierType
      return Object.keys(identifiers).find(key => identifiers[key] == id)
    },
    getWorkupLabel (type) {
      const workupTypes = reaction_pb.ReactionWorkup.ReactionWorkupType
      return Object.keys(workupTypes).find(key => workupTypes[key] == type).toLowerCase().replaceAll("_"," ")
    },
  },
  async mounted() {
    this.reaction = await this.getReactionData()
    this.reactionSummary = await this.getReactionSummary()
    this.loading = false
    console.log('schema',reaction_pb)
  }
}
</script>

<template lang="pug">
.main-reaction-view
  .section.summary(v-if='reactionSummary')
    .display(v-html='reactionSummary')
  .section(v-if='reaction?.identifiersList?.length')
    .title Identifiers
    .identifiers
      template(v-for='identifier in reaction.identifiersList')
        .value {{getReactionType(identifier.type)}}
        .value {{identifier.value}}
        .value {{identifier.details}}
  .section(v-if='reaction?.inputsMap?.length')
    .title Inputs
    .tabs
      .tab(
        v-for='(input, idx) in reaction.inputsMap'
        @click='inputsIdx = idx'
        :class='inputsIdx == idx ? "selected" : ""'
      ) {{input[0]}}
    .input
      .title Details
      .details
        template(v-for='key in Object.keys(displayInputs)')
          .label {{key.replaceAll(/(?=[A-Z])/g, ' ')}}
          .value {{displayInputs[key]}}
      .title Components
      .details
        template(v-for='component in reaction.inputsMap[inputsIdx][1].componentsList')
          CompoundView(
            :component='component'
          )
  .section(v-if='reaction?.setup')
    .title Setup
    .tabs
      template(
        v-for='tab in setupTabs'
      )
        .tab.capitalize(
          @click='setupTab = tab'
          :class='setupTab === tab ? "selected" : ""'
          v-if='tab !== "automation" || reaction.setup.is_automated'
        ) {{tab}}
    .details
      SetupView(
        :setup='reaction.setup'
        :display='setupTab'
      )
  .section(v-if='reaction?.conditions')
    .title Conditions
    .tabs
      template(
        v-for='tab in conditionTabs'
      )
        .tab.capitalize(
          @click='conditionTab = tab'
          :class='conditionTab === tab ? "selected" : ""'
          v-if='reaction.conditions[tab] || (tab === "other" && displayConditionsOther)'
        ) {{tab}}
    .details
      ConditionsView(
        :conditions='reaction.conditions'
        :display='conditionTab'
      )
  .section(v-if='reaction.notes')
    .title Notes
    .details
      NotesView(
        :notes='reaction.notes'
      )
  // TODO flesh out observations section
  .section(v-if='reaction.observationsList?.length')
    .title Observations
    .details
      ObservationsView(
        :observations='reaction.observationsList'
      )
  .section(v-if='reaction.workupsList?.length')
    .title Workups
    .tabs
      .tab.capitalize(
        v-for='(workup, idx) in reaction.workupsList'
        @click='workupsTab = idx'
        :class='workupsTab === idx ? "selected" : ""'
      ) {{getWorkupLabel(workup.type)}}
    .details
      WorkupsView(
        :workups='reaction.workupsList[workupsTab]'
      )

</template>

<style lang="sass" scoped>
.main-reaction-view
  margin-bottom: 2rem
.section
  width: calc(90vw)
  background-color: white
  border-radius: 0.25rem
  margin: 1rem auto 0
  padding: 1rem
  &.summary
    overflow-x: scroll
  .title
    font-weight: 700
    font-size: 1.5rem
    margin-bottom: 0.5rem
  .tabs
    display: flex
    column-gap: 0.5rem
    row-gap: 0.5rem
    margin-bottom: 0.5rem
    flex-wrap: wrap
    .tab
      padding: 0.5rem 1rem
      border-radius: 0.25rem
      border: 1px solid lightgrey
      cursor: pointer
      transition: 0.25s
      &.selected
        background-color: blue
        color: white
        border-color: blue
        cursor: default
      &.capitalize
        text-transform: capitalize
  .identifiers
    display: grid
    grid-template-columns: auto auto 1fr
    column-gap: 1rem
  .input
    .details
      display: grid
      grid-template-columns: auto 1fr
      column-gap: 1rem
      .label
        &:first-letter
          text-transform: uppercase

</style>