<script>
import EntityTable from '@/components/EntityTable'

export default {
  components: {
    EntityTable,
  },
  data() {
    return {
      loading: true,
      tableData: []
    }
  },
  mounted() {
    fetch("/api/fetch_datasets", {method: "GET"})
      .then(response => response.json())
      .then(data => {
        this.tableData = data
        this.loading = false
      })
  }
}
</script>

<template lang="pug">
#browse-main
  //- #datasets-table.table-container
  //-   .col.header Dataset ID
  //-   .col.header Name
  //-   .col.header Description
  //-   .col.header Size
  //-   template(
  //-     v-for='row in tableData'
  //-   )
  //-     .col {{row["Dataset ID"]}}
  //-     .col {{row.Name}}
  //-     .col.ellipses {{row.Description.length > 75 ? row.Description.substr(0,75)+"..." : row.Description}}
  //-     .col {{row.Size}}
  EntityTable(
    :tableData='tableData'
    title="",
    v-slot='{ entities }'
    :pagination='10'
    v-if='tableData.length'
  ) 
    .table-container
      .column.label Dataset ID
      .column.label Name
      .column.label Description
      .column.label Size
      template(
        v-for='(row, idx) in entities'
      )
        .column {{row["Dataset ID"]}}
        .column {{row.Name}}
        .column {{row.Description.length > 75 ? row.Description.substr(0,75)+"..." : row.Description}}
        .column {{row.Size}}
</template>

<style lang="sass" scoped>
@import '@/styles/table'
#browse-main
  padding: 1rem 0
  .table-container
    grid-template-columns: 1fr 1fr 1fr auto
</style>