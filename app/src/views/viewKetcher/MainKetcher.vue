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

<template lang="pug">
//- This is a simple view which emulates the html file and pulls in the js and css files for the ketcher interface
//- It expects a folder "ketcher" to be in the root of src with files extracted from:
//- https://github.com/epam/ketcher/releases/download/v2.5.1/ketcher-standalone-2.5.1.zip
#ketcher
  noscript You need to enable JavaScript to run this app.
  #root
</template>

<script>

export default {
  mounted() {
    // import ketcher js in mounted so the dom is ready
    import("@/ketcher/static/js/main.027562ee.js")
  }
}
</script>

<style lang="sass" scoped>
@import '@/ketcher/static/css/main.8e693d51.css'
</style>