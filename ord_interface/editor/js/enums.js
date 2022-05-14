/**
 * Copyright 2020 Open Reaction Database Project Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

goog.module('ord.enums');
goog.module.declareLegacyNamespace();

// Proto enums are referenced by reflection from "data-proto" HTML attributes.
/** @suppress {extraRequire} */
goog.require('proto.ord.ReactionRole.ReactionRoleType');
/** @suppress {extraRequire} */
goog.require('proto.ord.CompoundIdentifier.IdentifierType');
/** @suppress {extraRequire} */
goog.require('proto.ord.CompoundPreparation.PreparationType');
/** @suppress {extraRequire} */
goog.require('proto.ord.ElectrochemistryConditions.ElectrochemistryType');
/** @suppress {extraRequire} */
goog.require('proto.ord.FlowConditions.FlowType');
/** @suppress {extraRequire} */
goog.require('proto.ord.IlluminationConditions.IlluminationType');
/** @suppress {extraRequire} */
goog.require('proto.ord.PressureConditions.Atmosphere.AtmosphereType');
/** @suppress {extraRequire} */
goog.require('proto.ord.PressureConditions.Measurement.MeasurementType');
/** @suppress {extraRequire} */
goog.require('proto.ord.PressureConditions.PressureControl.PressureControlType');
/** @suppress {extraRequire} */
goog.require('proto.ord.ReactionIdentifier.IdentifierType');
/** @suppress {extraRequire} */
goog.require('proto.ord.ReactionInput.AdditionSpeed.AdditionSpeedType');
/** @suppress {extraRequire} */
goog.require('proto.ord.StirringConditions.StirringMethodType');
/** @suppress {extraRequire} */
goog.require('proto.ord.StirringConditions.StirringRate.StirringRateType');
/** @suppress {extraRequire} */
goog.require('proto.ord.TemperatureConditions.Measurement.MeasurementType');
/** @suppress {extraRequire} */
goog.require('proto.ord.TemperatureConditions.TemperatureControl.TemperatureControlType');
/** @suppress {extraRequire} */
goog.require('proto.ord.Vessel.VesselType');
/** @suppress {extraRequire} */
goog.require('proto.ord.VesselAttachment.VesselAttachmentType');
/** @suppress {extraRequire} */
goog.require('proto.ord.VesselMaterial.VesselMaterialType');
/** @suppress {extraRequire} */
goog.require('proto.ord.VesselPreparation.VesselPreparationType');
