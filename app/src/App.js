import './App.css';
import { BrowserRouter as Router, Routes, Route } from "react-router-dom"
import { Nav } from "./components/Nav"
import { Footer } from "./components/Footer"
// routes
import { Browse } from "./pages/Browse/index"

function App() {
  return (
    <Router>
      <div>
        <Nav />
        <Routes>
          <Route 
            exact 
            path="/"
          />
          <Route 
            exact 
            path="/browse" 
            element={<Browse/>}
          />
          <Route
            exact
            path="/search"
          />
          <Route 
            exact 
            path="/editor/datasets"
          />
        </Routes>
        <Footer />
      </div>
    </Router>
  );
}

export default App;
