# Copyright 2020 Open Reaction Database Project Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Tests for ord_interface.visualization.filters."""

from ord_schema.proto import reaction_pb2

from ord_interface.visualization import filters


def _pressure_with_setpoint() -> reaction_pb2.PressureConditions:
    pressure = reaction_pb2.PressureConditions()
    pressure.setpoint.value = 5
    pressure.setpoint.units = reaction_pb2.Pressure.ATMOSPHERE
    return pressure


def test_pressure_conditions_html_setpoint_without_atmosphere():
    # Regression: a setpoint with no atmosphere used to be dropped because the
    # setpoint block was gated on atmosphere.type rather than the setpoint.
    pressure = _pressure_with_setpoint()
    assert pressure.atmosphere.type == pressure.atmosphere.UNSPECIFIED
    assert "5 atm" in filters._pressure_conditions_html(pressure)


def test_pressure_conditions_html_setpoint_with_atmosphere():
    pressure = _pressure_with_setpoint()
    pressure.atmosphere.type = pressure.atmosphere.AIR
    result = filters._pressure_conditions_html(pressure)
    assert result == "in air (5 atm)"


def test_pressure_conditions_html_atmosphere_without_setpoint():
    pressure = reaction_pb2.PressureConditions()
    pressure.atmosphere.type = pressure.atmosphere.NITROGEN
    assert filters._pressure_conditions_html(pressure) == "under nitrogen"
