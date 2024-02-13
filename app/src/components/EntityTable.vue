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
export default {
  props: {
    tableData: Array,
    title: String,
    displaySearch: {type: Boolean, default: true}
  },
  watch: {
    pagination() {
      // set current page to 1 when pagination value changes to currentPage > lastPage
      this.currentPage = 1
    }
  },
  data() {
    return {
      entities: null,
      localFilters: [],
      search: '',
      currentPage: 1,
      loading: true,
      pagination: 10,
      searchString: "",
    }
  },
  computed: {
    filteredEntities() {
      if (!this.entities) return []
      let entities = this.entities

      if (this.searchString) {
        entities = entities.filter(entity => {
          // create a bool for each search param
          let matching = Array.apply(null, Array(this.searchArray.length))
          // check matching for each search param
          this.searchArray.forEach((param, pIdx) => {
            Object.keys(entity).forEach(key => {
              if (entity[key].toString().toLowerCase().includes(param.toLowerCase())) matching[pIdx] = true
            })
          })
          console.log('matching',matching.includes(undefined))
          return !matching.includes(undefined)
        })
      }

      return entities
    },
    paginatedEntities() {
      return this.filteredEntities.slice(this.pagiBottom, this.pagiTop)
    },
    pagiBottom() {
      return (this.currentPage - 1) * this.pagination
    },
    pagiTop() {
      return this.currentPage * this.pagination || 1
    },
    lastPage() {
      if (this.pagination && this.filteredEntities)
        return Math.ceil((this.filteredEntities.length || 1) / this.pagination)
      else return 1
    },
    pagiPrev() {
      return this.currentPage === 1 ? this.lastPage : this.currentPage - 1
    },
    pagiNext() {
      return this.currentPage === this.lastPage ? 1 : this.currentPage + 1
    },
    searchArray() {
      return this.searchString.split(" ")
    }
  },
  methods: {
  },
  async mounted() {
    this.entities = this.tableData
    // this.entities = Array.apply(null, Array(500)).map((val,idx) => {
    //   return {
    //     "Dataset ID": `${this.tableData[0]["Dataset ID"]}${Math.floor(Math.random() * idx)}`,
    //     Name: `${this.tableData[0].Name} ${Math.floor(Math.random() * idx)}`,
    //     Description: this.tableData[0].Description,
    //     Size: this.tableData[0].Size + Math.floor(Math.random() * idx),
    //   }
    // })
  }
}
</script>

<template lang="pug">
.table-main
  .table-container
    .header
      .title-holder
        .title {{ title }}
    .search-area(v-if='displaySearch')
      label(for="search") Search: 
      input(
        type="text"
        id="search"
        v-model='searchString'
      )
    .content
      .entities-holder
        slot(:entities='paginatedEntities')
      .pagination(
        v-if='pagination && entities'
      )
        .select Showing 
          select(
            name="pagination"
            id="pagination"
            v-model='pagination'
          )
            option(value=10) 10
            option(value=25) 25
            option(value=50) 50
            option(value=100) 100
          |  of {{filteredEntities.length}} entries.
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
@import '@/styles/vars'

.table-main
  width: 100%
  .table-container
    width: 100%
    .header
      .title-holder
        .title
          font-size: 2rem
          font-weight: 700
          margin-bottom: 1rem
    .search-area
      padding: 1rem
      margin: 0 5%
      input
        width: 30%
        min-width: 300px
    .content
      width: 100%
      height: fit-content
      .entities-holder
        // min-height: 14rem
        height: 100%
        display: block
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
        grid-template-columns: 1fr repeat(7, auto)
        // justify-items: end
        height: 100%
        align-items: center
        column-gap: 2rem
        padding: 1rem
        margin: 0 5%
        color: $linkblue
        .select
          color: black
        .paginav
          font-weight: 600
          display: inline-grid
          grid-column-gap: 0.33rem
          align-items: center
          cursor: pointer
          transition: .25s
          &:hover
            color: $hoverblue
          &.prev, &.next
            grid-template-columns: auto 1fr
          &.first, &.last
            grid-template-columns: auto auto 1fr
          &.no-click
            cursor: auto
          .chevron
            height: 0.75rem
        @media (max-width: 1000px)
          grid-template-columns: repeat(7, auto)
          row-gap: 1rem
          column-gap: 1rem
          .select
            grid-column: 1 / span 7
</style>