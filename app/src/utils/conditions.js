/**
 * Copyright 2023 Open Reaction Database Project Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import reaction_pb from "ord-schema"

export default {
  tempType(tempControlType) {
    const controlTypes = reaction_pb.TemperatureConditions.TemperatureControl
                             .TemperatureControlType
    return Object.keys(controlTypes)
        .find(key => controlTypes[key] == tempControlType)
  },
  tempSetPoint(setpoint) {
    if (!setpoint)
      return "None"
      const tempUnits = reaction_pb.Temperature.TemperatureUnit
      const tempUnit =
          Object.keys(tempUnits).find(key => tempUnits[key] == setpoint.units)
      return `${setpoint.value} (± ${setpoint.precision}) °${
          tempUnit.charAt(0)}`
  },
  pressureType(pressControlType) {
    const controlTypes =
        reaction_pb.PressureConditions.PressureControl.PressureControlType
    return Object.keys(controlTypes)
        .find(key => controlTypes[key] == pressControlType)
  },
  pressureSetPoint(setpoint) {
    if (!setpoint)
      return "None"
      const pressureUnits = reaction_pb.Pressure.PressureUnit
      const pressureUnit =
          Object.keys(pressureUnits)
              .find(key => pressureUnits[key] == setpoint.units)
      return `${setpoint.value} (± ${setpoint.precision}) ${
          pressureUnit.toLowerCase()}`
  },
  pressureAtmo(atmo) {
    const atmoTypes = reaction_pb.PressureConditions.Atmosphere.AtmosphereType
    const atmoType =
        Object.keys(atmoTypes).find(key => atmoTypes[key] == atmo.type)
    return `${atmoType}${atmo.details ? `, ${atmo.details}` : ""}`
  },
  stirType(stirType) {
    const stirTypes = reaction_pb.StirringConditions.StirringMethodType
    return Object.keys(stirTypes).find(key => stirTypes[key] == stirType)
  },
  stirRate(stirRate) {
    const stirRates = reaction_pb.StirringConditions.StirringRate
    return Object.keys(stirRates).find(key => stirRates[key] == stirRate)
  },
  illumType(illum) {
    const illumTypes = reaction_pb.IlluminationConditions.IlluminationType
    const illumType =
        Object.keys(illumTypes).find(key => illumTypes[key] == illum.type)
    return `${illumType}${illum.details ? `: ${illum.details}` : ""}`
  }
}
