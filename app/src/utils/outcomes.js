import reaction_pb from "ord-schema"

export default {
  formattedTime (timeData) {
    const timeUnits = reaction_pb.Time.TimeUnit
    const type = Object.keys(timeUnits).find(key => timeUnits[key] == timeData.units)
    return `${timeData.value} ${type.toLowerCase()}${timeData.units !== 0 ? "(s)" : ""}`
  }
}