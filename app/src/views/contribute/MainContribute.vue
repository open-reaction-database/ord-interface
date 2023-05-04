<script>

export default {
  components: {
  },
  data() {
    return {
      loading: true,
      tabs: [
        "Get Started",
        "Create",
        "Upload",
        "Enumerate"
      ],
      activeTab: "Get Started",
    }
  },
  computed: {
  },
  methods: {
  },
  async mounted() {
    this.loading = false
  },
}
</script>

<template lang="pug">
.main-contribute-view
  transition(name="fade")
    .loading(v-if='loading')
      LoadingSpinner
  transition(name="fade")
    .contribute-transition(v-if='!loading')
      .tabs
        .tab(
          v-for='tab in tabs'
          @click='activeTab = tab'
          :class='activeTab == tab ? "selected" : ""'
        ) {{tab}}
      transition(name="fade")
        .get-started(v-if='activeTab == "Get Started"')
          h3 Get Started
      transition(name="fade")
        .get-started(v-if='activeTab == "Create"')
          h3 Create
      transition(name="fade")
        .get-started(v-if='activeTab == "Upload"')
          h3 Upload
      transition(name="fade")
        .get-started(v-if='activeTab == "Enumerate"')
          h3 Enumerate
</template>

<style lang="sass" scoped>
@import "@/styles/vars"
@import "@/styles/tabs"
@import '@/styles/transition.sass'
</style>