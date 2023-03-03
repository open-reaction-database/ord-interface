<script>
import { reaction_pb } from "ord-schema"
import CompoundView from "./CompoundView"
import SetupView from "./SetupView"
import ConditionsView from "./ConditionsView"
import NotesView from "./NotesView"
import ObservationsView from "./ObservationsView"
import WorkupsView from "./WorkupsView"
import OutcomesView from "./OutcomesView"
import ProvenanceView from "./ProvenanceView"
import EventsView from "./EventsView"
import FloatingModal from "../../components/FloatingModal"

export default {
  components: {
    CompoundView,
    SetupView,
    ConditionsView,
    NotesView,
    ObservationsView,
    WorkupsView,
    OutcomesView,
    ProvenanceView,
    EventsView,
    FloatingModal,
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
      outcomesTab: 0,
      showRawReaction: false,
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
    },
    events() {
      const eventArray = []
      if (!this.reaction?.provenance?.recordCreated) return eventArray
      // add events to array
      eventArray.push(this.reaction.provenance.recordCreated)
      eventArray[0].details = "(record created)"
      eventArray.push(...this.reaction.provenance.recordModifiedList)
      // sort by date to be safe
      eventArray.sort((a,b) => {
        const dateA = new Date(a.time.value)
        const dateB = new Date(b.time.value)
        return dateA - dateB
      })
      return eventArray
    },
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
  template(v-if='reaction?.identifiersList?.length')
    .title Identifiers
    .section
      .identifiers
        template(v-for='identifier in reaction.identifiersList')
          .value {{getReactionType(identifier.type)}}
          .value {{identifier.value}}
          .value {{identifier.details}}
  template(v-if='reaction?.inputsMap?.length')
    .title Inputs
    .section
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
        .compound
          template(v-for='component in reaction.inputsMap[inputsIdx][1].componentsList')
            CompoundView(
              :component='component'
            )
  template(v-if='reaction?.setup')
    .title Setup
    .section
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
  template(v-if='reaction?.conditions')
    .title Conditions
    .section
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
  template(v-if='reaction.notes')
    .title Notes
    .section
      .details
        NotesView(
          :notes='reaction.notes'
        )
  // TODO flesh out observations section
  template(v-if='reaction.observationsList?.length')
    .title Observations
    .section
      .details
        ObservationsView(
          :observations='reaction.observationsList'
        )
  // TODO flesh out workups section
  template(v-if='reaction.workupsList?.length')
    .title Workups
    .section
      .tabs
        .tab.capitalize(
          v-for='(workup, idx) in reaction.workupsList'
          @click='workupsTab = idx'
          :class='workupsTab === idx ? "selected" : ""'
        ) {{getWorkupLabel(workup.type)}}
      .details
        WorkupsView(
          :workup='reaction.workupsList[workupsTab]'
        )
  template(v-if='reaction.outcomesList?.length')
    .title Outcomes
    .section
      .tabs
        .tab.capitalize(
          v-for='(outcome, idx) in reaction.outcomesList'
          @click='outcomesTab = idx'
          :class='outcomesTab === idx ? "selected" : ""'
        ) Outcome {{idx + 1}}
      .details
        OutcomesView(
          :outcome='reaction.outcomesList[outcomesTab]'
        )
  template(v-if='reaction.provenance')
    .title Provenance
    .section
      ProvenanceView(:provenance='reaction.provenance')
  template(v-if='events?.length')
    .title Record Events
    .section
      EventsView(:events='events')
  template(v-if='reaction')
    .title Full Record
    .section
      .full-record.button(@click='showRawReaction=true') View Full Record
    FloatingModal(
      v-if='showRawReaction'
      title="Raw Data"
      @closeModal='showRawReaction=false'
    )
      .data
        pre {{reaction}}
</template>

<style lang="sass" scoped>
@import "../../styles/tabs"
.main-reaction-view
  margin: 2rem 0
.title
  font-weight: 700
  font-size: 2rem
  margin-bottom: 0.5rem
.section, .title
  width: calc(95vw)
  margin: 0 auto
.section
  background-color: white
  border-radius: 0.25rem
  margin-bottom: 1rem
  padding: 1rem
  box-sizing: border-box
  .title
    font-size: 1.5rem
  &.summary
    overflow-x: scroll
  .identifiers
    display: grid
    grid-template-columns: auto auto 1fr
    column-gap: 1rem
  .input
    .details
      display: grid
      grid-template-columns: auto 1fr
      column-gap: 1rem
      margin-bottom: 1rem
      .label
        &:first-letter
          text-transform: uppercase
    .compound
      width: fit-content
  .full-record.button
    padding: 0.5rem 1rem
    background-color: blue
    border-radius: 0.25rem
    width: fit-content
    color: white
    cursor: pointer

</style>