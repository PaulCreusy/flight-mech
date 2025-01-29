"""
Module providing simple atmospheric models.
"""

###########
# Imports #
###########

# Dependencies #

import numpy as np
from scipy.optimize import minimize, Bounds

###########
# Classes #
###########

class LinearAtmosphere:
    """
    A very basic linear model for the atmosphere allowing to compute the density only.
    """

    rho_0 = 1.225  # kg.m-3

    @staticmethod
    def compute_sigma_from_altitude(z: float):
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

        sigma = (20 - z / 1e3) / (20 + z / 1e3)

        return sigma

    @staticmethod
    def compute_density_from_altitude(z: float):
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

        sigma = LinearAtmosphere.compute_sigma_from_altitude(z)
        rho = LinearAtmosphere.rho_0 * sigma

        return rho

    @staticmethod
    def compute_altitude_from_sigma(sigma: float):
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

        z = 20e3 * (1 - sigma) / (1 + sigma)

        return z


class StandardAtmosphere:
    """
    Standard atmosphere model to compute pressure, temperature and density with altitude.
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

        temperature = (288.15 + z * (216.65 - 288.15) / 11000) * (z < 11000) + \
            216.65 * (z >= 11000) * (z < 20000) + \
            (216.65 + (z - 20000) * (228.65 - 216.65) / 12000) * (z >= 20000) * (z < 32000) + \
            (228.65 + (z - 32000) * (270.65 - 228.65) / 12000) * \
            (z >= 32000) * (z < 47000)

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

        pressure = (101325 * np.power(1 - 22.588e-6 * z, 5.256)) * (z < 11000) + \
            (22632 * np.exp(-157.77e-6 * (z - 11000))) * (z >= 11000) * (z < 20000) + \
            (5474.9 * np.power(1 + 4.615e-6 * (z - 20000), 34.163)) * (z >= 20000) * (z < 32000) + \
            (868.014 * np.power(1 + 12.245e-6 * (z - 32000), -12.2)) * \
            (z >= 32000) * (z < 47000)

        return pressure

    @staticmethod
    def compute_density_from_altitude(z: float):
        """
        Compute the air density at a given altitude.

        Parameters
        ----------
        z : float
            Altitude in meters.

        Returns
        -------
        float
            Air density in kg.m-3.
        """

        # Compute pressure and temperature
        pressure = StandardAtmosphere.compute_pressure_from_altitude(z)
        temperature = StandardAtmosphere.compute_temperature_from_altitude(z)

        # Use perfect gas law to find density
        density = pressure / (StandardAtmosphere.r * temperature)

        return density

    @staticmethod
    def compute_sigma_from_altitude(z: float):
        """
        Compute the air density ratio at a given altitude.

        Parameters
        ----------
        z : float
            Altitude in meters.

        Returns
        -------
        float
            Air density ratio.
        """

        density = StandardAtmosphere.compute_density_from_altitude(z)
        sigma = density / StandardAtmosphere.rho_0

        return sigma

    @staticmethod
    def compute_sound_speed_from_altitude(z: float):
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

        temperature = StandardAtmosphere.compute_temperature_from_altitude(z)
        sound_speed = np.sqrt(StandardAtmosphere.gamma *
                              StandardAtmosphere.r * temperature)

        return sound_speed

    @staticmethod
    def compute_dynamic_viscosity_from_altitude(z: float):
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
        temperature = StandardAtmosphere.compute_temperature_from_altitude(z)

        # Compute mu
        mu = StandardAtmosphere.mu_0 * 0.083 * \
            np.power(temperature, 3 / 2) / (temperature + 110.4)

        return mu

    @staticmethod
    def compute_kinematic_viscosity_from_altitude(z: float):
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
        mu = StandardAtmosphere.compute_dynamic_viscosity_from_altitude(z)
        rho = StandardAtmosphere.compute_density_from_altitude(z)

        # Compute nu
        nu = mu / rho

        return nu

    @staticmethod
    def compute_altitude_from_sigma(sigma: float):
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
            return np.abs(StandardAtmosphere.compute_sigma_from_altitude(z) - sigma)

        # Find solution with scipy
        altitude = minimize(cost_function, 0, bounds=Bounds(0, 20e3)).x[0]

        return altitude
