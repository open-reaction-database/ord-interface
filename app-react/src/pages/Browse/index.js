import { useEffect, useState } from "react"
import "../../styles/browse.css"

const Browse = () => {
  const [tableData, setTableData] = useState([])

  useEffect(() => {
    fetch("client/api/fetch_datasets", {method: "GET"})
      .then(response => response.json())
      .then(data => {
        setTableData(data)
      })
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
          <tbody>
            {tableData.map((row) => (
              <tr>
                <td className="dataset_table_cell dataset_id">
                  {row["Dataset ID"]}
                </td>
                <td className="dataset_table_cell">
                  {row.Name}
                </td>
                <td className="dataset_table_cell ellipsis">
                  {row.Description}
                </td>
                <td className="dataset_table_cell">
                  {row.Size}
                </td>
              </tr>
            ))}
          </tbody>
      </table>
      <div className="table-container">
        <div className="row header">
          <div className="col header">Dataset ID</div>
          <div className="col header">Name</div>
          <div className="col header">Description</div>
          <div className="col header">Size</div>
        </div>
        {tableData.map((row) => (
          <div className="row">
            <div className="col">{row["Dataset ID"]}</div>
            <div className="col">{row.Name}</div>
            <div className="col">{row.Description}</div>
            <div className="col">{row.Size}</div>
          </div>
        ))}
      </div>
    </div>
  )
}

export default Browse