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
import DatasetCard from '@/components/DatasetCard'

export default {
  name: 'HomeView',
  components: {
    DatasetCard
  }
}
</script>

<template lang="pug">
.home
  .home-img
    img( alt="Ord logo" src="https://raw.githubusercontent.com/Open-Reaction-Database/ord-schema/main/logos/logo.svg")
  .dataset-showcase
    DatasetCard(
      :rank='0'
    )
    DatasetCard(
      :rank='1'
    )
    DatasetCard(
      :rank='2'
    )
    DatasetCard(
      :rank='3'
    )
</template>

<style lang="sass" scoped>
.home
  .home-img
    text-align: center
    width: 100%
    margin: 1% 20% 1% 20%
    max-width: 800px
    img
      width: 50%
  .dataset-showcase
    display: flex
    flex-direction: row
</style>
