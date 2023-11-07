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
  name: 'HomeView',
  components: {
  }
}
</script>

<template lang="pug">
.home
  .home-img
    img( alt="Ord logo" src="https://raw.githubusercontent.com/Open-Reaction-Database/ord-schema/main/logos/logo.svg")
</template>

<style lang="sass" scoped>
.home
  margin-top: 20vh
  .home-img
    width: 90%
    margin: auto
    max-width: 800px
    img
      width: 100%
</style>
