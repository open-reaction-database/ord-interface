<script>
import LoadingSpinner from '@/components/LoadingSpinner'

export default {
  components: {
    LoadingSpinner,
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
      activeTab: "Create",
      enumerateFiles: {
        template: {
          name: null,
          loading: false,
          value: null,
        },
        spreadsheet: {
          name: null,
          loading: false,
          value: null,
        }
      }
    }
  },
  computed: {
  },
  methods: {
    async setFile(e, fileType) {
      // converts uploaded file into useable string
      const files = e.target.files || e.dataTransfer.files
      if (!files.length) return console.error('No file')
      this.enumerateFiles[fileType].loading = true
      this.enumerateFiles[fileType].name = files[0].name
      const fileReader = new FileReader()
      fileReader.onload = readerEvent => {
        this.enumerateFiles[fileType].value = readerEvent.target.result
        this.enumerateFiles[fileType].loading = false
      }
      if (fileType == "template")
        fileReader.readAsText(files[0])
      else if (fileType == "spreadsheet")
        fileReader.readAsDataURL(files[0])
    },
    submitEnumerate() {
      if (this.enumerateFiles.template.loading || this.enumerateFiles.spreadsheet.loading)
        return alert("Files are still processing. Please try again in a moment.")
      else if (!this.enumerateFiles.template.value || !this.enumerateFiles.spreadsheet.value)
        return alert("You must upload a file for the template and spreadsheet before submitting.")
      // send enumerated files to api for upload
      const xhr = new XMLHttpRequest();
      xhr.open('POST', '/editor-api/dataset/enumerate');
      xhr.onload = function () {
        if (xhr.status === 200) {
          location.reload();
        } else {
          alert(`Error: ${xhr.response}`)
          console.error(xhr.response)
        }
      }
      xhr.send(JSON.stringify({
        'template_string': this.enumerateFiles.template.value,
        'spreadsheet_data': this.enumerateFiles.spreadsheet.value,
        'spreadsheet_name': this.enumerateFiles.spreadsheet.name
      }));
    }
  },
  async mounted() {
    // if this is user's first time, default page is get started
    const notFirstTime = localStorage.getItem('notFirstTime')
    if (!notFirstTime) {
      localStorage.setItem("notFirstTime", "false")
      this.activeTab = "Get Started"
    }
    this.loading = false
  },
}
</script>

<template lang="pug">
.main-contribute
  transition(name="fade")
    .loading(v-if='loading')
      LoadingSpinner
  transition(name="fade")
    .contribute-menu(v-if='!loading')
      .title {{activeTab}}
      .contribute-container
        .tabs
          .tab(
            v-for='tab in tabs'
            @click='activeTab = tab'
            :class='activeTab == tab ? "selected" : ""'
          ) {{tab}}
        transition(name="fade")
          .get-started(v-if='activeTab == "Get Started"')
            .copy Check out these 
              a(href="https://www.youtube.com/playlist?list=PLyoEVAlMb276aRRa4xLNRAzbMPRlNb7VI") tutorial videos 
              | for a guide on how to use the ORD Editor.
            .tutorial-videos
              .video
                .subtitle Using the interactive web editor (updated 1/2021)
                iframe(
                  width='560'
                  height='315'
                  src='https://www.youtube.com/embed/5J-2j8aBXMo'
                  title='YouTube video player'
                  frameborder='0'
                  allow='accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share'
                  allowfullscreen
                )
              .video
                .subtitle Submitting your first Dataset
                iframe(
                  width='560'
                  height='315'
                  src='https://www.youtube.com/embed/sAdSKKdO9Gs'
                  title='YouTube video player'
                  frameborder='0'
                  allow='accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share'
                  allowfullscreen
                )
            .copy Please send any questions, comments, or issues to 
              a(
                href='https://mail.google.com/mail/?view=cm&fs=1&tf=1&to=help@open-reaction-database.org'
                target="_blank"
              ) help@open-reaction-database.org 
              | or create a new issue on the 
              a(
                href='https://github.com/Open-Reaction-Database/ord-interface/issues'
                target="_blank"
              ) ord-interface GitHub repository
              | .
        transition(name="fade")
          .get-started(v-if='activeTab == "Create"')
        transition(name="fade")
          .get-started(v-if='activeTab == "Upload"')
        transition(name="fade")
          .enumerate(v-if='activeTab == "Enumerate"')
            .subtitle Upload files for enumeration:
            .file-picker
              .input
                label(for='template') Template filename:
                input#template(
                  type='file'
                  accept='.pbtxt'
                  v-on:change='(e) => setFile(e,"template")'
                )
              .input
                label(for='spreadsheet') Spreadsheet filename:
                input#spreadsheet(
                  type='file'
                  accept='.csv,.xls,.xlsx'
                  v-on:change='(e) => setFile(e,"spreadsheet")'
                )
            .submit
              button#enumerate_submit(@click='submitEnumerate') Submit Enumeration Upload
            .copy You can use a Reaction template to enumerate a Dataset based on a spreadsheet of values. For more information, see the 
              a(
                href='https://docs.open-reaction-database.org/en/latest/guides/templates.html'
                target="_blank"
              ) documentation
              | .
            .copy
              b NOTE:
              |  Large dataset enumerations (thousands of reactions) may result in a browser timeout. If this happens, please send an email to 
              a(href='mailto:help@open-reaction-database.org') help@open-reaction-database.org
              |  and attach your template and spreadsheet files. Alternatively, you may use the 
              a(href='https://github.com/Open-Reaction-Database/ord-schema/blob/main/ord_schema/scripts/enumerate_dataset.py') programmatic interface
              |  to enumerate the dataset locally.

</template>

<style lang="sass" scoped>
@import "@/styles/vars"
@import "@/styles/tabs"
@import '@/styles/transition.sass'
.main-contribute
  margin: 0 5%
  padding: 1rem
  .title
    font-size: 2rem
    font-weight: 700
    margin-bottom: 1rem
  .subtitle
    font-size: 1.5rem
  .contribute-container
    background-color: white
    border-radius: 0.25rem
    padding: 1rem
  .get-started
    .copy
      margin: 3rem 0 2rem
    .tutorial-videos
      display: flex
      flex-wrap: wrap
      column-gap: 1rem
      row-gap: 1rem
  .enumerate
    padding-top: 1rem
    .file-picker
      padding: 1rem 0
      display: flex
      flex-wrap: wrap
      row-gap: 1rem
      .input
        label
          margin-right: 0.5rem
    .submit
      margin-bottom: 2rem
    .copy
      margin-top: 1rem
</style>