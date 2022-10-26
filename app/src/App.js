import './App.css';
import { BrowserRouter as Router, Routes, Route, Link } from "react-router-dom"

function App() {
  return (
    <Router>
      <div>
        <Routes>
          <Route exact path="/">

          </Route>
          <Route exact path="/client/browse">

          </Route>
          <Route exact path="/search">

          </Route>
          <Route exact path="/editor/datasets">

          </Route>
        </Routes>
      </div>
    </Router>
  );
}

export default App;
