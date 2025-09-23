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

export const amountObj =
    (amount) => {
      if (!amount)
        return {}
      let units = {};
      let unitCategory = "";
      // determine unit type
      if (amount.moles) {
        units = reaction_pb.Moles.MolesUnit
        unitCategory = "moles"
      } else if (amount.volume) {
        units = reaction_pb.Volume.VolumeUnit
        unitCategory = "volume"
      } else if (amount.mass) {
        units = reaction_pb.Mass.MassUnit
        unitCategory = "mass"
      } else if (amount.unmeasured) {
        return { unitAmount: "", unitCategory: "unmeasured" }
      }
      const unitVal = amount[unitCategory].units
      const amountVal = amount[unitCategory].value
      const precision = amount[unitCategory].precision
      return {
        unitAmount: amountVal,
            unitType: Object.keys(units).find(key => units[key] == unitVal),
            unitCategory: unitCategory,
      }
    }

export const amountStr = (amountObj) => {
  // takes an amountObj from above amountObj function
  if (!amountObj.unitAmount || !amountObj.unitType)
    return ""
    return `${Math.round(amountObj.unitAmount * 1000) / 1000} ${
        amountObj.unitType.toLowerCase()}`
}
