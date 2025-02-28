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
}
</script>

<template lang="pug">
nav.navbar.navbar-expand-lg.bg-light
  .container
    a.navbar-brand(href="/")
      img(
        src="https://raw.githubusercontent.com/Open-Reaction-Database/ord-schema/main/logos/logo.svg"
        alt="ORD Logo"
        height="30"
      )
    #navbarNav.collapse.navbar-collapse
      .navbar-nav
        .nav-item
          router-link.nav-link(:to='{name: "browse"}') Browse
          router-link.nav-link(:to='{name: "search"}') Search
        .nav-item
          a.nav-link(href="https://open-reaction-database.org/editor/") Contribute
          a.nav-link(href="https://docs.open-reaction-database.org") Docs
        .nav-item
          router-link.nav-link(:to='{name: "about"}') About
</template>

<style lang="sass" scoped>
@import @/styles/vars

nav
  background-color: white !important
  .container
    display: flex
    flex-wrap: inherit
    align-items: center
    width: 90%
    margin: auto
    padding: 0.5rem 0
    .navbar-brand
      padding: 0.5rem 0
    #navbarNav
      flex-direction: row
      .navbar-nav
        margin-left: 1rem
        .nav-item
          display: flex
          .nav-link
            padding: 0.5rem
            text-decoration: none
            transition: .15s
            color: $linkblue
            &:hover
              color: $hoverblue
</style>
