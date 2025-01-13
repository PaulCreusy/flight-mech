"""
Module providing a simple atmospheric model.
"""

#############
# Constants #
#############

rho_0 = 1.225  # kg.m-3

#############
# Functions #
#############

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

def compute_air_density_from_altitude(z: float):
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

    sigma = compute_sigma_from_altitude(z)
    rho = rho_0 * sigma

    return rho

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
