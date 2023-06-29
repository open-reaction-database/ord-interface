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
      formattedResults: []
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
        return `Yield: ${yieldObj.percentage.value}%`
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
    template(
      v-for='row in entities'
    )
      .reaction-link
        .row
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
              .yield {{getYield(row.data.outcomesList[0].productsList[0].measurementsList)}}
              .conditions {{conditionsAndDuration(row.data).join("; ")}}
              .smile(v-if='row.data.outcomesList[0].productsList[0].identifiersList.length')
                .value {{productIdentifier(row.data.outcomesList[0].productsList[0].identifiersList[0])}}
                CopyButton(
                  :textToCopy='row.data.outcomesList[0].productsList[0].identifiersList[0].value'
                )
            .col
              .creator Uploaded by {{row.data.provenance.recordCreated.person.name}}, {{row.data.provenance.recordCreated.person.organization}}
              .date Uploaded on {{new Date(row.data.provenance.recordCreated.time.value).toLocaleDateString()}}
              .publication 
                a(
                  :href='row.data.provenance.publicationUrl'
                  target="_blank"
                ) Publication URL

</template>

<style lang="sass" scoped>
@import '@/styles/vars.sass'
.search-results-main
  .action-button-holder
    margin: -2.5rem 0 1rem // bring button row inline with title without having to pass too much into EntityTable
    display: flex
    justify-content: flex-end
    button
  .reaction-link
    text-decoration: none
    .row
      background-color: white
      border-radius: 0.25rem
      padding: 1rem
      margin-bottom: 1rem
      transition: 0.25s
      &:hover
        box-shadow: 0 0 5px $darkgrey
        cursor: pointer
      .reaction-table
        color: black
        overflow-x: wrap
      .info
        display: grid
        grid-template-columns: repeat(2, 1fr)
        row-gap: 1rem
        column-gap: 1rem
        margin-top: 1rem
        .col.full
          grid-column: 1/3
          button
            font-size: 1.2rem
  @media (max-width: 1000px)
    margin-top: 2.5rem
</style>