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

      // ensure that there are filters and at least one is engaged
      if (
        this.localFilters &&
        this.localFilters.length &&
        this.localFilters.find((lf) => lf.engaged)
      ) {
        entities = this.localFilters.reduce(
          (acc, filter) => (filter.engaged ? acc.filter(filter.func) : acc),
          entities
        )
      }

      // ensure that there are dropdowns and at least one has chosen values
      // if (
      //   this.localDropdowns &&
      //   this.localDropdowns.length &&
      //   this.localDropdowns.find((ld) => ld.chosen.length)
      // ) {
      //   entities = this.localDropdowns.reduce((acc, ld) => {
      //     if (!ld.chosen.length) return acc

      //     return acc.filter((entity) => ld.func(entity, ld.chosen))
      //   }, entities)
      // }

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
    readyFilters() {
      if (!(this.filters && this.filters.length)) return

      this.localFilters = this.filters.map((filter, index) => ({
        toggleExclusive: () => {
          // If selecting this filter should deselect another filter
          if (filter.exclusive.length)
            filter.exclusive.forEach((partner) => {
              this.localFilters[partner].engaged = false
            })
          this.localFilters[index].engaged = true
        },
        ...filter,
      }))
    },
    handleFilterClick(filter) {
      if (filter.exclusive.length) filter.toggleExclusive()
      else filter.engaged = !filter.engaged
    },
  },
  async mounted() {
    this.entities = this.tableData
    console.log('entities',this.entities)
    // await this.readyFilters()
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
    .filters(v-if='localFilters && localFilters.length')
      .filter(
        v-for='filter in localFilters',
        @click='handleFilterClick(filter)',
        :class='[{ engaged: filter.engaged }, filter.class && filter.class()]'
      ) {{ filter.copy(filter.engaged) }}
    .content
      .entities-holder
        slot(:entities='filteredEntities')
      .pagination(
        v-if='pagination && entities && entities.length > pagination'
      )
        .prev.paginav(@click='currentPage = pagiPrev') 
          img.chevron(src='/img/arrowL.png')
          span.word Previous
        .next.paginav(@click='currentPage = pagiNext') 
          span.word Next
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
    .filters
      margin-bottom: 1rem
      display: flex
      user-select: none
      flex-wrap: wrap
      .filter
        cursor: pointer
        margin: 0 1rem .5rem 0
        // border: thin solid $yellow
        padding: .5rem 1rem
        // border-radius: .25rem
        width: fit-content
        font-size: 12px
        box-shadow: 0 0 .5rem 0 #a0a0a0
        &.underline
          box-shadow: none
          &:active
            box-shadow: none
        &.engaged
          background-color: grey
          color: white
          &.underline
            background-color: transparent
            color: black
            border-bottom: 5px solid grey
        &:active
          box-shadow: 0 0 .5rem 0 rgba(0,0,0,.4) inset
    .search-area
      padding: 1rem 0
    // .dropdowns
    //   position: relative
    //   display: grid
    //   grid-template-columns: 1fr 1fr
    //   column-gap: 1rem
    //   margin-bottom: 1rem
    //   .dropdown-holder
    //     display: grid
    //     align-items: center
    //     column-gap: .25rem
    //     .helper-buttons
    //       .button
    //         color: white
    //         // background: $mild-blue
    //         padding: .25rem .5rem
    //         // border-radius: .25rem
    //         cursor: pointer
    //         font-size: 12px

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
        grid-template-columns: 1fr auto
        justify-items: end
        height: 100%
        align-items: end
        column-gap: 2rem
        padding: 1rem
        .paginav
          font-weight: 600
          display: inline-grid
          grid-template-columns: auto 1fr
          grid-column-gap: 0.33rem
          align-items: center
          &:hover
            cursor: pointer
          .chevron
            height: 0.75rem
</style>