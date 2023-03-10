

<script>
export default {
}
</script>

<template lang="pug">
nav.navbar.navbar-expand-lg.bg-light
  .container
    a.navbar-brand(href="https://open-reaction-database.org")
      img(
        src="https://raw.githubusercontent.com/Open-Reaction-Database/ord-schema/main/logos/logo.svg" 
        alt="ORD Logo"
        height="30"
      )
    //- button(
    //-   className="navbar-toggler" 
    //-   type="button" 
    //-   data-bs-toggle="collapse" 
    //-   data-bs-target="#navbarNav"
    //-   aria-controls="navbarNav" 
    //-   aria-expanded="false"
    //-   aria-label="Toggle navigation"
    //- )
    //-   span.navbar-toggler-icon
    #navbarNav.collapse.navbar-collapse
      .navbar-nav
        .nav-item
          router-link.nav-link(to="/browse") Browse
          router-link.nav-link(to="/search") Search
        //- li.nav-item
        //-   <Link className="nav-link" to="/editor">Contribute</Link>
        .nav-item
          a.nav-link(href="https://docs.open-reaction-database.org") Docs
</template>

<style lang="sass" scoped>
@import @/styles/vars

nav
  background-color: white
  .container
    display: flex
    flex-wrap: inherit
    align-items: center
    // justify-content: space-between
    width: 90%
    margin: auto
    padding: 0.5rem 0
    .navbar-brand
      padding: 0.5rem 0
    #navbarNav
      flex-direction: row
      .navbar-nav
        display: flex
        margin-left: 1rem
        .nav-item
          .nav-link
            padding: 0.5rem
            text-decoration: none
            transition: .15s
            color: $linkblue
            &:hover
              color: $hoverblue
</style>