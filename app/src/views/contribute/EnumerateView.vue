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
  data() {
    return {
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
}
</script>

<template lang="pug">
.enumerate
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
    a(
      href='mailto:help@open-reaction-database.org'
      target="_blank"
    ) help@open-reaction-database.org
    |  and attach your template and spreadsheet files. Alternatively, you may use the 
    a(
      href='https://github.com/Open-Reaction-Database/ord-schema/blob/main/ord_schema/scripts/enumerate_dataset.py'
      target="_blank"
    ) programmatic interface
    |  to enumerate the dataset locally.

</template>

<style lang="sass" scoped>
.subtitle
  font-size: 1.5rem
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