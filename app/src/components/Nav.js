import { Link } from "react-router-dom"

const Nav = () => {
  return (
    <nav class="navbar navbar-expand-lg bg-light">
      <div class="container">
        <a class="navbar-brand" href="https://open-reaction-database.org">
          <img 
            src="https://raw.githubusercontent.com/Open-Reaction-Database/ord-schema/main/logos/logo.svg" 
            alt=""
            height="30px" />
        </a>
        <button 
          class="navbar-toggler" 
          type="button" 
          data-bs-toggle="collapse" 
          data-bs-target="#navbarNav"
          aria-controls="navbarNav" 
          aria-expanded="false"
          aria-label="Toggle navigation"
        >
          <span class="navbar-toggler-icon"></span>
        </button>
        <div class="collapse navbar-collapse" id="navbarNav">
          <ul class="navbar-nav">
            <li class="nav-item">
              <Link class="nav-link active" aria-current="page" to="/browse">Browse</Link>
            </li>
            <li class="nav-item">
              <Link class="nav-link" to="/search">Search</Link>
            </li>
            <li class="nav-item">
              <Link class="nav-link" to="/editor">Contribute</Link>
            </li>
            <li class="nav-item">
              <Link class="nav-link" to="https://docs.open-reaction-database.org">Docs</Link>
            </li>
          </ul>
        </div>
      </div>
    </nav>
  )
}

export default Nav