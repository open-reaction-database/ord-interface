<script>
export default {
  props: {
    tableData: Array,
    pagination: Number,
  },
  data() {
    return {
      entities: null,
      localFilters: [],
      search: '',
      currentPage: 1,
      loading: true,
    }
  },
  computed: {
    filteredEntities() {if (!this.entities) return []
      let entities = this.entities

      // if (this.search) {
      //   entities = entities.filter(entity => isBeingSearched(this.search, (entity[this.searchOptions.field] || '')))
      // }

      if (this.pagination > 0) {
        entities = entities.slice(this.pagiBottom, this.pagiTop)
      }

      return entities
    },
    pagiBottom() {
      return (this.currentPage - 1) * this.pagination
    },
    pagiTop() {
      return this.currentPage * this.pagination || 1
    },
    lastPage() {
      if (this.pagination && this.entities)
        return Math.ceil((this.entities.length || 1) / this.pagination)
      else return 1
    },
    pagiPrev() {
      return this.currentPage === 1 ? this.lastPage : this.currentPage - 1
    },
    pagiNext() {
      return this.currentPage === this.lastPage ? 1 : this.currentPage + 1
    },
  },
  methods: {
  },
  async mounted() {
    // this.entities = this.tableData
    this.entities = Array.apply(null, Array(100)).map(idx => this.tableData[0])
  }
}
</script>

<template lang="pug">
.table-main
  .table-container
    .header
      .title-holder
        .title.dotdotdot {{ title }}
        .subtitle.dotodotdot {{ subtitle }}
    .search-area(v-if='searchOptions && searchOptions.field')
      .input-holder hello world
        //- pro-input(
        //-   v-model='search',
        //-   :options='{ title: `Search by ${searchOptions.by}` }'
        //- )
    .content
      .entities-holder
        slot(:entities='filteredEntities')
      .pagination(
        v-if='pagination && entities && entities.length > pagination'
      )
        .paginav.first(@click='currentPage = 1')
          img.chevron(src='/img/arrowL.png')
          img.chevron(src='/img/arrowL.png')
          span.word First
        .prev.paginav(@click='currentPage = pagiPrev') 
          img.chevron(src='/img/arrowL.png')
          span.word Previous
        .paginav(
          @click='currentPage = currentPage > 1 ? currentPage - 1 : currentPage'
          :class='currentPage > 1 ? "" : "no-click"'
        )
          span.word {{ currentPage > 1 ? currentPage - 1 : "..."}}
        .paginav.no-click
          span.word {{ currentPage }}
        .paginav(
          @click='currentPage = currentPage < lastPage ? currentPage + 1 : currentPage'
          :class='currentPage < lastPage ? "" : "no-click"'
        )
          span.word {{ currentPage < lastPage ? currentPage + 1 : "..."}}
        .next.paginav(@click='currentPage = pagiNext') 
          span.word Next
          img.chevron(src='/img/arrowR.png')
        .paginav.last(@click='currentPage = lastPage')
          span.word Last
          img.chevron(src='/img/arrowR.png')
          img.chevron(src='/img/arrowR.png')
</template>

<style lang="sass" scoped>
.table-main
  width: 100%
  .table-container
    width: 100%
    .header
      display: grid
      grid-template-rows: auto auto
      grid-template-columns: auto auto
      row-gap: 1rem
      .breadcrumbs
        grid-column: 1/3
        grid-row: 1/2
      .title-holder
        display: grid
        grid-template-columns: 100%
        grid-template-rows: auto auto
        .subtitle
          font-weight: bold
    .search-area
      padding: 1rem 0

    .content
      width: 100%
      height: fit-content
      .entities-holder
        min-height: 14rem
        height: 100%
        // overflow-y: auto
        // overflow-x: hidden
        display: block
        // box-shadow: 0 0 .25rem 0 $lightgrey inset
        width: 100%
        > *
          padding-left: 1rem
          padding-right: 1rem
        .entity.no-interact
          display: grid
          cursor: default
          grid-template-columns: 1fr
          column-gap: 1rem
          align-items: center
          padding: 1rem
          &:hover
            background-color: transparent
      .pagination
        display: grid
        grid-template-columns: 1fr repeat(6, auto)
        justify-items: end
        height: 100%
        align-items: end
        column-gap: 2rem
        padding: 1rem
        margin-right: 5%
        .paginav
          font-weight: 600
          display: inline-grid
          grid-column-gap: 0.33rem
          align-items: center
          cursor: pointer
          &.prev, &.next
            grid-template-columns: auto 1fr
          &.first, &.last
            grid-template-columns: auto auto 1fr
          &.no-click
            cursor: auto
          .chevron
            height: 0.75rem
</style>