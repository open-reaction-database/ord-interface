import { useEffect } from 'react'

const Browse = () => {

  useEffect(() => {
    const tableData = [];
    console.log('hello')
    fetch("http://localhost:5000/")
      .then(response => response.json())
      .then(data => console.log('data',data))
    return 
  }, [])

  return (
    <div id="datasets" className="container mt-3">
      <table id="datasets_table" className="display table">
          <thead>
          <tr>
              <th>Dataset ID</th>
              <th>Name</th>
              <th>Description</th>
              <th>Size</th>
          </tr>
          </thead>
      </table>
    </div>
  )
}

export default Browse