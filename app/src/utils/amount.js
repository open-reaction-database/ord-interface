import { reaction_pb } from "ord-schema"

export const amountObj = (amount) => {
  if (!amount) return {}
  let units = {}
  let unitCategory = ""
  // determine unit type
  if (amount.moles) {
    units = reaction_pb.Moles.MolesUnit
    unitCategory="moles"
  } else if (amount.volume) {
    units = reaction_pb.Volume.VolumeUnit
    unitCategory = "volume"
  } else if (amount.mass) {
    units = reaction_pb.Mass.MassUnit
    unitCategory = "mass"
  } else if (amount.unmeasured) {
    return {unitAmount: "", unitCategory: "unmeasured"}
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
  if (!amountObj.unitAmount || !amountObj.unitType) return ""
  return `${Math.round(amountObj.unitAmount * 1000) / 1000} ${amountObj.unitType.toLowerCase()}`
}