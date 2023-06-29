import reaction_pb from "ord-schema"

export default {
  tempType (tempControlType) {
    const controlTypes = reaction_pb.TemperatureConditions.TemperatureControl.TemperatureControlType
    return Object.keys(controlTypes).find(key => controlTypes[key] == tempControlType)
  },
  tempSetPoint (setpoint) {
    // const setpoint = this.conditions.temperature.setpoint
    if (!setpoint) return "None"
    const tempUnits = reaction_pb.Temperature.TemperatureUnit
    const tempUnit = Object.keys(tempUnits).find(key => tempUnits[key] == setpoint.units)
    const tempSymbols = {
      CELSIUS: "C",

    }
    return `${setpoint.value} (± ${setpoint.precision}) °${tempUnit.charAt(0)}`
  },
  pressureType (pressControlType) {
    const controlTypes = reaction_pb.PressureConditions.PressureControl.PressureControlType
    return Object.keys(controlTypes).find(key => controlTypes[key] == pressControlType)
  },
  pressureSetPoint (setpoint) {
    if (!setpoint) return "None"
    const pressureUnits = reaction_pb.Pressure.PressureUnit
    const pressureUnit = Object.keys(pressureUnits).find(key => pressureUnits[key] == setpoint.units)
    return `${setpoint.value} (± ${setpoint.precision}) ${pressureUnit.toLowerCase()}`
  },
  pressureAtmo (atmo) {
    const atmoTypes = reaction_pb.PressureConditions.Atmosphere.AtmosphereType
    const atmoType = Object.keys(atmoTypes).find(key => atmoTypes[key] == atmo.type)
    return `${atmoType}${atmo.details ? `, ${atmo.details}` : ""}`
  },
  stirType (stirType) {
    const stirTypes = reaction_pb.StirringConditions.StirringMethodType
    return Object.keys(stirTypes).find(key => stirTypes[key] == stirType)
  },
  illumType (illum) {
    const illumTypes = reaction_pb.IlluminationConditions.IlluminationType
    const illumType = Object.keys(illumTypes).find(key => illumTypes[key] == illum.type)
    return `${illumType}${illum.details ? `: ${illum.details}` : ""}`
  }
}
