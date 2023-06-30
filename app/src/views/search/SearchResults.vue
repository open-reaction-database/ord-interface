<script>
import EntityTable from '@/components/EntityTable'
import LoadingSpinner from '@/components/LoadingSpinner'
import conditionUtil from '@/utils/conditions'
import outcomesUtil from '@/utils/outcomes'
import reaction_pb from 'ord-schema'
import CopyButton from '@/components/CopyButton'

export default {
  props: {
    searchResults: Array,
  },
  components: {
    EntityTable,
    LoadingSpinner,
    CopyButton
  },
  data() {
    return {
      formattedResults: [],
      selectedReactions: [],
    }
  },
  methods: {
    getReactionTables() {
      this.formattedResults.forEach(result => {
        fetch(`/api/render/${result.reaction_id}`)
          .then(response => response.json())
          .then(responseData => {
            result.reactionTable = responseData
          })
      })
    },
    downloadResults() {
      // create .pb download of search results
      const xhr = new XMLHttpRequest();
      xhr.open('POST', 'api/download_results');
      xhr.responseType = "blob";
      xhr.onload = () => {
        if (xhr.status === 200) {
          const url = URL.createObjectURL(xhr.response);
          const link = document.createElement('a');
          link.href = url;
          link.download = "ord_search_results.pb.gz"
          link.click();
          // https://stackoverflow.com/a/56547307.
          setTimeout(() => {
            URL.revokeObjectURL(url);
            link.remove();
          }, 100);
        }
      };
      xhr.setRequestHeader('Content-Type', 'application/json');

      // format request json to expected key/value in api
      const requestJson = this.formattedResults.map(result => {return {"Reaction ID": result.reaction_id}})
      xhr.send(JSON.stringify(requestJson));
    },
    getYield(measurements) {
      const yieldObj = measurements.find(m => m.type == 3) // ord-schema type 3 == "YIELD"
      if (yieldObj.percentage) {
        return `${yieldObj.percentage.value}%`
      } else {
        return ""
      }
    },
    conditionsAndDuration(reaction) {
      const details = []
      // get temp
      const tempSetPoint = conditionUtil.tempSetPoint(reaction.conditions.temperature.setpoint)
      if (tempSetPoint !== "None")
        details.push(`at ${tempSetPoint}`)

      // get Pressure
      const pressureSetPoint = conditionUtil.pressureSetPoint(reaction.conditions.pressure.setpoint)
      if (pressureSetPoint !== "None")
        details.push(`under ${pressureSetPoint}`)

      // get duration
      const formattedTime = outcomesUtil.formattedTime(reaction.outcomesList[0].reactionTime)
      details.push(`for ${formattedTime}`)

      return details
    },
    productIdentifier(identifier) {
      const identifierTypes = reaction_pb.CompoundIdentifier.CompoundIdentifierType
      const identifierType = Object.keys(identifierTypes).find(key => identifierTypes[key] == identifier.type)
      return `${identifierType}: ${identifier.value}`
    },
    updateSelectedReactions(event) {
      if (event.target.checked) {
        this.selectedReactions.push(event.target.value)
      } else {
        let idx = this.selectedReactions.indexOf(event.target.value)
        if (idx !== -1) {
          this.selectedReactions.splice(idx, 1)
        }
      }
    },
    goToViewSelected() {
      this.$router.push({ name: 'selected-set', query: {reaction_ids: this.selectedReactions}})
    }
  },
  async mounted() {
    this.formattedResults = this.searchResults
    this.getReactionTables()
  },
}
</script>

<template lang="pug">
.search-results-main
  EntityTable(
    :tableData='formattedResults'
    title="Search Results",
    v-slot='{ entities }'
    v-if='formattedResults.length'
    :displaySearch='false'
  ) 
    .action-button-holder
      button(
        :disabled='!formattedResults.length'
        @click='downloadResults'
      ) Download Results
    .reaction-container(
      v-for='row in entities'
    )
      .row(:class='selectedReactions.includes(row.reaction_id) ? "selected" : ""')
        .select
          input(
            type="checkbox"
            :id='"select_"+row.reaction_id'
            :value='row.reaction_id'
            :checked='selectedReactions.includes(row.reaction_id)'
            @change='updateSelectedReactions($event)'
          )
          label(:for='"select_"+row.reaction_id') Select reaction
        .reaction-table(
          v-html='row.reactionTable'
          v-if='row.reactionTable'
        )
        LoadingSpinner(v-else)
        .info
          .col.full
            router-link(
              :to='{ name: "reaction-view", params: {reactionId: row.reaction_id}}'
            ) 
              button View Full Details
          .col
            .yield Yield: {{getYield(row.data.outcomesList[0].productsList[0].measurementsList)}}
            .conditions Conditions: {{conditionsAndDuration(row.data).join("; ")}}
            .smile(v-if='row.data.outcomesList[0].productsList[0].identifiersList.length')
              CopyButton(
                :textToCopy='row.data.outcomesList[0].productsList[0].identifiersList[0].value'
              )
              .value {{productIdentifier(row.data.outcomesList[0].productsList[0].identifiersList[0])}}
          .col
            .creator Uploaded by {{row.data.provenance.recordCreated.person.name}}, {{row.data.provenance.recordCreated.person.organization}}
            .date Uploaded on {{new Date(row.data.provenance.recordCreated.time.value).toLocaleDateString()}}
            .doi DOI: {{row.data.provenance.doi}}
            .publication 
              a(
                :href='row.data.provenance.publicationUrl'
                target="_blank"
              ) Publication URL
  transition(
    name='fade'
  )
    .view-selected-container(v-if='selectedReactions.length')
      .view-selected-button(@click='goToViewSelected') View {{selectedReactions.length}} selected reactions
</template>

<style lang="sass" scoped>
@import '@/styles/vars.sass'
@import '@/styles/transition.sass'
.search-results-main
  .action-button-holder
    margin: -2.5rem 0 1rem // bring button row inline with title without having to pass too much into EntityTable
    display: flex
    justify-content: flex-end
    button
  .reaction-container
    text-decoration: none
    .row
      background-color: white
      border-radius: 0.25rem
      padding: 0.75rem
      margin-bottom: 1rem
      transition: 0.25s
      border: 4px solid white
      &:hover
        box-shadow: 0 0 5px $darkgrey
      .select
        input, label
          cursor: pointer
      .reaction-table
        color: black
        overflow-x: wrap
      .info
        display: grid
        grid-template-columns: repeat(2, 50%)
        row-gap: 0.5rem
        column-gap: 1rem
        margin-top: 1rem
        .col.full
          grid-column: 1/3
          button
            font-size: 1.2rem
        .col
          *
            margin-top: 0.5rem
        .smile
          display: flex
          column-gap: 0.5rem
          align-items: center
          margin-top: 0.25rem
          .value
            width: 100%
            white-space: nowrap
            overflow: hidden
            text-overflow: ellipsis
      &.selected
        border-color: $linkblue
  .view-selected-container
    position: fixed
    bottom: 2rem
    right: 6.5rem
    .view-selected-button
      padding: 0.5rem 1rem
      color: white
      background-color: $linkblue
      border-radius: 0.25rem
      cursor: pointer
  @media (max-width: 1000px)
    margin-top: 2.5rem

</style>