"""
Module providing simple atmospheric models.
"""

###########
# Imports #
###########

# Python imports #

from typing import Literal
from abc import ABC, abstractmethod

# Dependencies #

import numpy as np
from scipy.optimize import minimize, Bounds

# Local imports #

from flight_mech._common import plot_graph

#############
# Constants #
#############

VARIABLE_TO_UNIT = {
    "pressure": "Pa",
    "temperature": "K",
    "density": "kg.m-3",
    "sound_speed": "m.s-1",
    "dynamic_viscosity": "kg.m-1.s-1",
    "kinematic_viscosity": "m2.s-1"
}

#############
# Functions #
#############

def compute_air_sound_speed(temperature: float):
    """
    Compute the air sound speed at the given temperature, assuming it is a perfect gas.

    Parameters
    ----------
    temperature : float
        Temperature of the air in K.

    Returns
    -------
    float
        Sound speed in m.s-1.
    """

    sound_speed = np.sqrt(StandardAtmosphere.gamma *
                          StandardAtmosphere.r * temperature)

    return sound_speed

###########
# Classes #
###########

class AtmosphereModel(ABC):
    """
    Atmosphere abstract class.
    """

    rho_0: float

    @staticmethod
    @abstractmethod
    def compute_sigma_from_altitude(z: float) -> float:
        """
        Compute the density coefficient sigma.

        Parameters
        ----------
        z : float
            Altitude in meters.

        Returns
        -------
        float
            Density coefficient.
        """
        pass

    @classmethod
    @abstractmethod
    def compute_density_from_altitude(self, z: float) -> float:
        """
        Compute the air density.

        Parameters
        ----------
        z : float
            Altitude in meters.

        Returns
        -------
        float
            Air density in kg.m-3
        """
        pass

    @staticmethod
    @abstractmethod
    def compute_altitude_from_sigma(sigma: float) -> float:
        """
        Compute the altitude corresponding to the given density coefficient.

        Parameters
        ----------
        sigma : float
            Density coefficient.

        Returns
        -------
        float
            Altitude in meters.
        """
        pass

    def plot_graph(
            self,
            variable: Literal["pressure", "temperature", "density", "sound_speed", "dynamic_viscosity"],
            min_altitude: float = 0.,
            max_altitude: float = 20000.,
            nb_points: int = 100,
            **kwargs):
        """
        Plot the graph of evolution of the given variable in the atmosphere.

        Parameters
        ----------
        variable : Literal["pressure", "temperature", "density", "sound_speed", "dynamic_viscosity"]
            Name of the variable to plot.
        min_altitude : float, optional
            Minimum altitude on the plot in meters, by default 0.
        max_altitude : float, optional
            Maximum altitude on the plot in meters, by default 20000.
        nb_points : int, optional
            Number of points in the plot.

        Note
        ----
        For more details on the optional arguments, please check flight_mech._common.plot_graph.
        """

        # Define an altitude array
        altitude_array = np.linspace(min_altitude, max_altitude, nb_points)

        # Compute the variable array
        variable_array = getattr(
            self, f"compute_{variable}_from_altitude")(altitude_array)

        # Plot
        plot_graph(
            x_array=variable_array,
            y_array=altitude_array,
            title=f"{variable.capitalize()} graph",
            y_label="Altitude in meters",
            x_label=f"{variable.capitalize()} [{VARIABLE_TO_UNIT[variable]}]",
            **kwargs
        )

class ConstantAtmosphere(AtmosphereModel):
    """
    A constant atmosphere model. Used for test purposes only.
    """

    rho_0 = 1.225  # kg.m-3

    @staticmethod
    def compute_sigma_from_altitude(z: float):
        return 1.

    @staticmethod
    def compute_altitude_from_sigma(sigma: float):
        return 0.

    @classmethod
    def compute_density_from_altitude(self, z: float):
        return self.rho_0


class LinearAtmosphere(AtmosphereModel):
    """
    A very basic linear model for the atmosphere allowing to compute the density only.
    """

    rho_0 = 1.225  # kg.m-3

    @staticmethod
    def compute_sigma_from_altitude(z: float):
        sigma = (20 - z / 1e3) / (20 + z / 1e3)
        return sigma

    @classmethod
    def compute_density_from_altitude(self, z: float):
        sigma = self.compute_sigma_from_altitude(z)
        rho = self.rho_0 * sigma
        return rho

    @staticmethod
    def compute_altitude_from_sigma(sigma: float):
        z = 20e3 * (1 - sigma) / (1 + sigma)
        return z


class StandardAtmosphere(AtmosphereModel):
    """
    International Standard Atmosphere model to compute pressure, temperature and density with altitude.

    Notes
    -----
    For more details, see: https://en.wikipedia.org/wiki/International_Standard_Atmosphere.
    For the original paper, see: https://www.digitaldutch.com/atmoscalc/US_Standard_Atmosphere_1976.pdf.
    """

    gamma = 1.4
    r = 287.058  # J.kg-1.K-1
    rho_0 = 1.225  # kg.m-3
    mu_0 = 17.26e-6  # kg.m-1.s-1

    @staticmethod
    def compute_temperature_from_altitude(z: float):
        """
        Compute the temperature at a given altitude.

        Parameters
        ----------
        z : float
            Altitude in meters.

        Returns
        -------
        float
            Air temperature in K.
        """

        temperature = (288.15 + z * -6.5 / 1000) * (z < 11000) + \
            216.65 * (z >= 11000) * (z < 20000) + \
            (216.65 + (z - 20000) * 1. / 1000) * (z >= 20000) * (z < 32000) + \
            (228.65 + (z - 32000) * 2.8 / 1000) * \
            (z >= 32000) * (z < 47000) + \
            270.15 * (z >= 47000) * (z < 51000) + \
            (270.15 + (z - 51000) * -2.8 / 1000) * (z >= 51000) * (z < 71000) + \
            (214.15 + (z - 71000) * -2. / 1000) * (z >= 71000) * (z <= 86000)

        return temperature

    @staticmethod
    def compute_pressure_from_altitude(z: float):
        """
        Compute the pressure at a given altitude.

        Parameters
        ----------
        z : float
            Altitude in meters.

        Returns
        -------
        float
            Air pressure in Pa.
        """

        pressure = (z < 11000) * (101325 * np.power(np.abs(1 - 22.588e-6 * z), 5.256)) + \
            (z >= 11000) * (z < 20000) * (22632 * np.exp(-157.77e-6 * (z - 11000))) + \
            (5474.9 * np.power(1 - 4.615e-6 * (z - 20000), 34.163)) * (z >= 20000) * (z < 32000) + \
            (868.014 * np.power(1 + 12.245e-6 * (z - 32000), -12.2)) * \
            (z >= 32000) * (z < 47000) + \
            (110.906 * np.exp(-126.293e-6 * (z - 47000))) * (z >= 47000) * (z < 51000) + \
            (66.939 * np.power((270.65) / (270.65 - 2.8 * (z - 51000) / 1000), -12.2)) * \
            (z >= 51000) * (z < 71000) + \
            (3.9564 * np.power((214.65) / (214.65 - 2. * (z - 71000) / 1000), -17.09)) * \
            (z >= 71000) * (z <= 86000)

        return pressure

    @classmethod
    def compute_density_from_altitude(self, z: float):

        # Compute pressure and temperature
        pressure = self.compute_pressure_from_altitude(z)
        temperature = self.compute_temperature_from_altitude(z)

        # Use perfect gas law to find density
        density = pressure / (self.r * temperature)

        return density

    @classmethod
    def compute_sigma_from_altitude(self, z: float):

        density = self.compute_density_from_altitude(z)
        sigma = density / self.rho_0

        return sigma

    @classmethod
    def compute_sound_speed_from_altitude(self, z: float):
        """
        Compute the sound speed at a given altitude.

        Parameters
        ----------
        z : float
            Altitude in meters.

        Returns
        -------
        float
            Sound speed in m.s-1.
        """

        temperature = self.compute_temperature_from_altitude(z)
        sound_speed = compute_air_sound_speed(temperature)

        return sound_speed

    @classmethod
    def compute_dynamic_viscosity_from_altitude(self, z: float):
        """
        Compute the dynamic viscosity at a given altitude.

        Parameters
        ----------
        z : float
            Altitude in meters.

        Returns
        -------
        float
            Dynamic viscosity in kg.m-1.s-1.
        """

        # Compute temperature
        temperature = self.compute_temperature_from_altitude(z)

        # Compute mu
        mu = self.mu_0 * 0.083 * \
            np.power(temperature, 3 / 2) / (temperature + 110.4)

        return mu

    @classmethod
    def compute_kinematic_viscosity_from_altitude(self, z: float):
        """
        Compute the kinematic viscosity at a given altitude.

        Parameters
        ----------
        z : float
            Altitude in meters.

        Returns
        -------
        float
            Kinematic viscosity in m2.s-1.
        """

        # Compute mu and rho
        mu = self.compute_dynamic_viscosity_from_altitude(z)
        rho = self.compute_density_from_altitude(z)

        # Compute nu
        nu = mu / rho

        return nu

    @classmethod
    def compute_altitude_from_sigma(self, sigma: float):
        """
        Compute the altitude from a given density ratio. This only works for altitude below 20km.

        Parameters
        ----------
        z : float
            Density ratio.

        Returns
        -------
        float
            Altitude in m.
        """

        # Define a cost function to minimize
        def cost_function(z):
            return np.abs(self.compute_sigma_from_altitude(z) - sigma)

        # Find solution with scipy
        altitude = minimize(cost_function, 0, bounds=Bounds(0, 20e3)).x[0]

        return altitude
