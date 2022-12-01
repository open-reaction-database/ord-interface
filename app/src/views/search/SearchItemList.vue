<script>
export default {
  props: {
    title: String,
    itemList: Array,
  },
  data () {
    return {
      mutatedList: [],
      itemToAdd: "",
    }
  },
  methods: {
    updateList() {
      this.mutatedList.push(this.itemToAdd)
      this.itemToAdd = ""
      this.$emit('update:itemList', this.mutatedList)
    },
    deleteItem(idx) {
      this.mutatedList.splice(idx,1)
      this.$emit('update:itemList', this.mutatedList)
    },
  },
  mounted() {
    this.mutatedList = this.itemList || []
  }
}
</script>

<template lang="pug">
.item-list
  .title {{title}}
  .list 
    template(v-for='(item, idx) in mutatedList')
      .text {{item}}
      .delete
        button(@click='deleteItem(idx)')
          i.material-icons delete
  .input
    input(
      type="text"
      v-model='itemToAdd'
    )
    button(@click='updateList')
      i.material-icons add
      span Add
</template>

<style lang="sass" scoped>
.item-list
  .title
    font-size: 1.25rem
    font-weight: 700
  .list
    display: grid
    grid-template-columns: auto 1fr
    column-gap: 0.5rem
    row-gap: 0.5rem
    margin-bottom: 1rem
  .input
    display: flex
  button
    display: flex
    align-items: center
    margin-left: 0.5rem
    i
      font-size: 1.1rem
</style>