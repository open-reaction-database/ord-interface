import { Link } from "react-router-dom"

const Nav = () => {
  return (
    <nav className="navbar navbar-expand-lg bg-light">
      <div className="container">
        <a className="navbar-brand" href="https://open-reaction-database.org">
          <img 
            src="https://raw.githubusercontent.com/Open-Reaction-Database/ord-schema/main/logos/logo.svg" 
            alt=""
            height="30px" />
        </a>
        <button 
          className="navbar-toggler" 
          type="button" 
          data-bs-toggle="collapse" 
          data-bs-target="#navbarNav"
          aria-controls="navbarNav" 
          aria-expanded="false"
          aria-label="Toggle navigation"
        >
          <span className="navbar-toggler-icon"></span>
        </button>
        <div className="collapse navbar-collapse" id="navbarNav">
          <ul className="navbar-nav">
            <li className="nav-item">
              <Link className="nav-link active" aria-current="page" to="/browse">Browse</Link>
            </li>
            <li className="nav-item">
              <Link className="nav-link" to="/search">Search</Link>
            </li>
            <li className="nav-item">
              <Link className="nav-link" to="/editor">Contribute</Link>
            </li>
            <li className="nav-item">
              <Link className="nav-link" to="https://docs.open-reaction-database.org">Docs</Link>
            </li>
          </ul>
        </div>
      </div>
    </nav>
  )
}

export default Nav