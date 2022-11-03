<template lang="pug">
#datasets.container.mt-3
  #datasets-table.table-container
    .row.header
      .col.header Dataset ID
      .col.header Name
      .col.header Description
      .col.header Size
    .row(
      v-for='row in tableData'
    )
      .col {{row["Dataset ID"]}}
      .col {{row.Name}}
      .col {{row.Description}}
      .col {{row.Size}}
</template>

<script>
export default {
  data() {
    return {
      tableData: []
    }
  },
  mounted() {
    fetch("/api/fetch_datasets", {method: "GET"})
      .then(response => response.json())
      .then(data => {
        this.tableData = data
      })
  }
}
</script>

<style lang="sass" scoped>

</style>