<template lang="pug">
#browse-main
  #datasets-table.table-container
    .col.header Dataset ID
    .col.header Name
    .col.header Description
    .col.header Size
    template(
      v-for='row in tableData'
    )
      .col {{row["Dataset ID"]}}
      .col {{row.Name}}
      .col.ellipses {{row.Description.length > 75 ? row.Description.substr(0,75)+"..." : row.Description}}
      .col {{row.Size}}
</template>

<script>
export default {
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

<style lang="sass" scoped>
@import '@/styles/table'
#browse-main
  padding: 1rem 0
  #datasets-table
    grid-template-columns: 1fr 1fr 1fr auto
</style>