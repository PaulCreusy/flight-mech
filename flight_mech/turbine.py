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
from flight_mech.fuel import FuelTypes
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

class TurbojetSingleBody:
    """
    Class to define a turbojet with single flux, single body with a constant Cp coefficient.
    """

    ambient_pressure = 101325  # Pa
    ambient_temperature = 285  # K
    T4_max = 1700  # K
    compressor_efficiency = 0.86
    turbine_efficiency = 0.9
    M0 = 0
    OPR_design = 10
    max_reference_surface_mass_flow_rate_4_star = 241.261  # kg.s-1.m-2
    A4_star = 1e-2  # m-2
    mode: Literal["design", "operation"] = "design"
    fuel = FuelTypes.KEROSENE

    @property
    def P0(self):
        P0 = self.ambient_pressure * \
            np.power((1 + (GAMMA - 1) / 2 * np.power(self.M0, 2)),
                     GAMMA / (GAMMA - 1))
        return P0

    @property
    def T0(self):
        T0 = self.ambient_temperature * \
            (1 + (GAMMA - 1) / 2 * np.power(self.M0, 2))
        return T0

    @property
    def P1(self):
        return self.P0

    @property
    def T1(self):
        return self.T0

    @property
    def P2(self):
        return self.P0

    @property
    def T2(self):
        return self.T0

    @property
    def P3(self):
        return self.P2 * self.OPR_design

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
        return self.T4_max

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
        M8 = np.sqrt(
            (np.power(self.P8 / self.Ps8, (GAMMA - 1) / GAMMA) - 1) * (2 / (GAMMA - 1)))
        return M8

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
        A8_star = self.W8R / self.max_reference_surface_mass_flow_rate_4_star
        return A8_star

    @property
    def A8(self):
        A8 = self.A8_star * (1 / self.M8) * np.power((2 / (GAMMA + 1)) * (
            1 + (GAMMA - 1) / 2 * np.power(self.M8, 2)), (GAMMA + 1) / (2 * (GAMMA - 1)))
        return A8

    @property
    def v8(self):
        v8 = self.M8 * compute_air_sound_speed(self.Ts8)
        return v8

    @property
    def thrust(self):
        thrust = self.v8 * self.W8
        return thrust

    @property
    def Wf(self):
        Wf = (1 / self.fuel.lower_heating_value) * \
            self.W4 * cp * (self.T4 - self.T3)
        return Wf

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

    def plot_graph(self, variable: Literal["pressure", "temperature", "mass_flow"], **kwargs):

        # Extract code of the variable
        code = VARIABLE_TO_CODE[variable]

        # Allocate list for the x and y values
        values_list = []
        id_list = []

        # Extract the values
        for i in range(1, 9):
            id_list.append(i)
            values_list.append(self.__getattribute__(code + str(i)))

        # Plot
        plot_graph(
            x_array=id_list,
            y_array=values_list,
            title=f"{variable.capitalize()} graph",
            x_label="Section id",
            y_label=f"{variable.capitalize()} [{VARIABLE_TO_UNIT[variable]}]",
            **kwargs
        )
