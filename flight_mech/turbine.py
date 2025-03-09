"""
Module to define turbo-reactors models and compute their characteristics.
"""

###########
# Imports #
###########

# Python imports #

from typing import Literal

# Dependencies #

import numpy as np
from scipy.optimize import brute

# Local imports #

from flight_mech.atmosphere import (
    StandardAtmosphere,
    compute_air_sound_speed
)
from flight_mech.fuel import FuelTable
from flight_mech._common import plot_graph

#############
# Constants #
#############

# Define reference quantities for computations
REFERENCE_PRESSURE = 101325  # Pa
REFERENCE_TEMPERATURE = 288.15  # K
GAMMA = StandardAtmosphere.gamma
cp = (GAMMA * StandardAtmosphere.r) / (GAMMA - 1)

VARIABLE_TO_CODE = {
    "pressure": "P",
    "temperature": "T",
    "mass_flow": "W"
}

VARIABLE_TO_UNIT = {
    "pressure": "Pa",
    "temperature": "K",
    "mass_flow": "kg.s-1"
}


###########
# Classes #
###########

class _Turbojet:
    """
    Base class to define turbojets.
    """

    M0: float = 0
    _ambient_pressure: float = 101325  # Pa
    _ambient_temperature: float = 285  # K
    _altitude: float | None = None
    _ambient_conditions_from_altitude: bool = False

    atmosphere_model = StandardAtmosphere

    mode: Literal["design", "operation"] = "design"
    fuel = FuelTable.KEROSENE

    @property
    def ambient_pressure(self):
        return self._ambient_pressure

    @ambient_pressure.setter
    def ambient_pressure(self, value):
        self._ambient_pressure = value
        self._ambient_conditions_from_altitude = False
        self._altitude = None

    @property
    def ambient_temperature(self):
        return self._ambient_temperature

    @ambient_temperature.setter
    def ambient_temperature(self, value):
        self._ambient_temperature = value
        self._ambient_conditions_from_altitude = False
        self._altitude = None

    @property
    def altitude(self):
        return self._altitude

    @altitude.setter
    def altitude(self, value: float):
        self._altitude = value
        self._ambient_temperature = self.atmosphere_model.compute_temperature_from_altitude(
            value)
        self._ambient_pressure = self.atmosphere_model.compute_pressure_from_altitude(
            value)
        self._ambient_conditions_from_altitude = True

    def plot_graph(self, variable: Literal["pressure", "temperature", "mass_flow"], **kwargs):
        """
        Plot the graph of evolution of the given variable in the turbojet.

        Parameters
        ----------
        variable : Literal["pressure", "temperature", "mass_flow"]
            Name of the variable to plot.

        Note
        ----
        For more details on the optional arguments, please check flight_mech._common.plot_graph.
        """

        # Extract code of the variable
        code = VARIABLE_TO_CODE[variable]

        # Allocate list for the x and y values
        values_list = []
        id_list = []

        # Extract the values
        for i in range(1, 9):
            current_code = code + str(i)
            if current_code in self.__dir__():
                id_list.append(i)
                values_list.append(self.__getattribute__(current_code))

        # Plot
        plot_graph(
            x_array=id_list,
            y_array=values_list,
            title=f"{variable.capitalize()} graph",
            x_label="Section id",
            y_label=f"{variable.capitalize()} [{VARIABLE_TO_UNIT[variable]}]",
            **kwargs
        )


class TurbojetSingleBody(_Turbojet):
    """
    Class to define a turbojet with single flux, single body with a constant Cp coefficient.
    """

    compressor_efficiency = 0.86
    turbine_efficiency = 0.9

    T4_max = 1700  # K
    OPR_design = 10
    max_reference_surface_mass_flow_rate_4_star = 241.261  # kg.s-1.m-2
    A4_star = 1e-2  # m-2

    # Define operation variables
    _T4_instruction = T4_max
    current_OPR = OPR_design

    @property
    def T4_instruction(self):
        return self._T4_instruction

    @T4_instruction.setter
    def T4_instruction(self, value):
        if value < 0:
            raise ValueError(
                f"T4 cannot accept negative values ({value}), please provide only positive values.")
        self._T4_instruction = value

    @property
    def P0(self):
        P0 = self.ambient_pressure * \
            np.power((1 + (GAMMA - 1) / 2 * np.power(self.M0, 2)),
                     GAMMA / (GAMMA - 1))
        return P0

    @property
    def V0(self):
        return self.M0 * compute_air_sound_speed(self.ambient_temperature)

    @property
    def T0(self):
        T0 = self.ambient_temperature * \
            (1 + (GAMMA - 1) / 2 * np.power(self.M0, 2))
        return T0

    @property
    def W0(self):
        return self.W3

    @property
    def P1(self):
        return self.P0

    @property
    def T1(self):
        return self.T0

    @property
    def W1(self):
        return self.W3

    @property
    def P2(self):
        return self.P0

    @property
    def T2(self):
        return self.T0

    @property
    def W2(self):
        return self.W3

    @property
    def P3(self):
        if self.mode == "design":
            return self.P2 * self.OPR_design
        elif self.mode == "operation":
            return self.P2 * self.current_OPR

    @property
    def T3(self):
        T3 = self.T2 * (1 + (1 / self.compressor_efficiency) *
                        (np.power(self.P3 / self.P2, (GAMMA - 1) / GAMMA) - 1))
        return T3

    @property
    def W3(self):
        W3 = self.W4 - self.Wf
        return W3

    @property
    def P4(self):
        P4 = 0.95 * self.P3
        return P4

    @property
    def T4(self):
        if self.mode == "design":
            return self.T4_max
        elif self.mode == "operation":
            return self.T4_instruction

    @property
    def W4R(self):
        W4R = self.max_reference_surface_mass_flow_rate_4_star * self.A4_star
        return W4R

    @property
    def W4(self):
        W4 = self.W4R * (self.P4 / REFERENCE_PRESSURE) / \
            np.sqrt(self.T4 / REFERENCE_TEMPERATURE)
        return W4

    @property
    def Wf(self):
        Wf = (1 / self.fuel.lower_heating_value) * \
            self.W4 * cp * (self.T4 - self.T3)
        return Wf

    @property
    def P5(self):
        P5 = self.P4 * np.power(1 - (1 / self.turbine_efficiency)
                                * (1 - (self.T5 / self.T4)), GAMMA / (GAMMA - 1))
        return P5

    @property
    def T5(self):
        T5 = self.T4 - (self.T3 - self.T2)
        return T5

    @property
    def W5(self):
        return self.W4

    @property
    def P8(self):
        return self.P5

    @property
    def T8(self):
        return self.T5

    @property
    def Ps8(self):
        return self.ambient_pressure

    @property
    def M8(self):
        if self.mode == "design":
            M8 = np.sqrt(
                (np.power(self.P8 / self.Ps8, (GAMMA - 1) / GAMMA) - 1) * (2 / (GAMMA - 1)))
            return M8
        elif self.mode == "operation":
            # Use this property because all the computations need to be performed in design mode
            return self._get_design_variable("M8")

    @property
    def Ts8(self):
        Ts8 = self.T8 * np.power(1 + (GAMMA - 1) / 2 *
                                 np.power(self.M8, 2), -1)
        return Ts8

    @property
    def W8(self):
        return self.W5

    @property
    def W8R(self):
        W8R = self.W8 * np.sqrt(self.T8 / REFERENCE_TEMPERATURE) / \
            (self.P8 / REFERENCE_PRESSURE)
        return W8R

    @property
    def A8_star(self):
        if self.mode == "design":
            A8_star = self.W8R / self.max_reference_surface_mass_flow_rate_4_star
            return A8_star
        elif self.mode == "operation":
            # Use this property because all the computations need to be performed in design mode
            return self._get_design_variable("A8_star")

    @property
    def A8(self):
        if self.mode == "design":
            A8 = self.A8_star * (1 / self.M8) * np.power((2 / (GAMMA + 1)) * (
                1 + (GAMMA - 1) / 2 * np.power(self.M8, 2)), (GAMMA + 1) / (2 * (GAMMA - 1)))
            return A8
        elif self.mode == "operation":
            # Use this property because all the computations need to be performed in design mode
            return self._get_design_variable("A8")

    @property
    def V8(self):
        v8 = self.M8 * compute_air_sound_speed(self.Ts8)
        return v8

    @property
    def thrust(self):
        self._check_temperatures_positivity()
        thrust = self.V8 * self.W8 + self.A8 * \
            (self.Ps8 - self.ambient_pressure)
        return thrust

    @property
    def fuel_consumption(self):
        return self.Wf

    @property
    def global_efficiency(self):
        global_efficiency = self.thrust * self.V0 / \
            (self.Wf * self.fuel.lower_heating_value)
        return global_efficiency

    @property
    def propulsive_efficiency(self):
        propulsive_efficiency = self.thrust * self.V0 / \
            (self.W0 * (np.square(self.V8) - np.square(self.V0)) / 2)
        return propulsive_efficiency

    @property
    def thermal_efficiency(self):
        thermal_efficiency = (self.W0 * (np.square(self.V8) - np.square(self.V0)) / 2) /\
            (self.Wf * self.fuel.lower_heating_value)
        return thermal_efficiency

    def _check_temperatures_positivity(self):
        for i in range(1, 9):
            current_code = f"T{i}"
            if current_code not in self.__dir__():
                continue
            current_value = self.__getattribute__(current_code)
            if current_value < 0:
                raise ValueError(
                    f"{current_code} is negative ({current_value}), the domain is outside its domain of validity.")

    def _get_design_variable(self, variable: str) -> float:
        # Raise error if already in design mode
        if self.mode == "design":
            raise ValueError(
                f"You are currently in design mode. Please use the property directly instead.")

        # Store the previous mode
        previous_mode = self.mode

        # Switch to design mode
        self.mode = "design"

        # Compute M8
        value = self.__getattribute__(variable)

        # Switch back to previous mode
        self.mode = previous_mode

        return value

    def tune_A4_star_for_desired_thrust(self, desired_thrust: float, min_A4_star: float = 1e-4, max_A4_star: float = 5e-1):
        """
        Tune the value of A4* to obtain the desired thrust in the operating conditions.

        Parameters
        ----------
        desired_thrust : float
            Desired thrust in N.
        min_A4_star : float, optional
            Minimal value for A4*, by default 1e-4
        max_A4_star : float, optional
            Maximum value for A4*, by default 5e-1
        """

        # Define a cost function
        def cost_function(A4_star):
            self.A4_star = A4_star
            return np.abs(desired_thrust - self.thrust)

        # Solve by brute force
        res = brute(cost_function, [(min_A4_star, max_A4_star)])

        # Update A4*
        self.A4_star = res[0]

    def tune_current_OPR(self):
        """
        Tune the current OPR value to the T4 instruction and operating conditions.

        Raises
        ------
        ValueError
            Raise error if not in operation mode.
        """

        # Raise error if not in operation mode
        if self.mode != "operation":
            raise ValueError(
                "The OPR convergence process can only be used in operating mode. In design mode, the current OPR is always considered to be the design OPR.")

        # Define a cost function
        def cost_function(current_OPR):
            self.current_OPR = current_OPR
            cost = np.abs((self.W8R / self.A8_star) -
                          (self.W4R / self.A4_star))
            return cost

        # Solve by brute force
        res = brute(cost_function, [(0, self.OPR_design)])

        # Update OPR
        self.current_OPR = res[0]
