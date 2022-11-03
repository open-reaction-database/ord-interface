<template lang="pug">
nav.navbar.navbar-expand-lg.bg-light
    .container
      a.navbar-brand(href="https://open-reaction-database.org")
        img(
          src="https://raw.githubusercontent.com/Open-Reaction-Database/ord-schema/main/logos/logo.svg" 
          alt=""
          height="30px"
        )
      button(
        className="navbar-toggler" 
        type="button" 
        data-bs-toggle="collapse" 
        data-bs-target="#navbarNav"
        aria-controls="navbarNav" 
        aria-expanded="false"
        aria-label="Toggle navigation"
      )
        span.navbar-toggler-icon
      #navbarNav.collapse.navbar-collapse
        ul.navbar-nav
          li.nav-item
            router-link.nav-link.active(aria-current="page" to="/browse") Browse
          //- li.nav-item
          //-   <Link className="nav-link" to="/search">Search</Link>
          //- li.nav-item
          //-   <Link className="nav-link" to="/editor">Contribute</Link>
          li.nav-item
            a.nav-link(href="https://docs.open-reaction-database.org") Docs
</template>

<script>
export default {
  
}
</script>

<style lang="sass" scoped>

</style>