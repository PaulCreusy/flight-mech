"""
Module to define functions to compute aerodynamic effects.
"""

###########
# Imports #
###########

# Python imports #

from typing import Literal

# Dependencies #

import numpy as np
from scipy.integrate import solve_ivp, simpson
import matplotlib.pyplot as plt

#############
# Constants #
#############

# Define the transition Reynolds number
turbulent_reynolds = 5e5

#############
# Functions #
#############

def check_turbulent_transition(x: float | np.ndarray, velocity: float | np.ndarray, nu: float):
    """
    Check if a velocity distribution along a plate has a turbulent transition.

    Parameters
    ----------
    x : float | np.ndarray
        Array of coordinates along the plate.
    velocity : float | np.ndarray
        Array of velocity values.
    nu : float
        Kinematic viscosity.

    Raises
    ------
    ValueError
        Raise error if turbulent transition.
    """

    # Compute the reynolds number
    reynolds = np.array(x * velocity / nu)

    # Check transition
    if (reynolds > turbulent_reynolds).any():
        raise ValueError(
            f"The boundary layer becomes turbulent. The Reynolds number at the end of the plate is {'{:3e}'.format(reynolds)}.")

def compute_blasius_linear_drag(x_array: np.ndarray, velocity_array: np.ndarray, rho: float, nu: float, return_array: bool = False):
    """
    Compute the linear drag under the flat rectangle, uniform velocity hypothesis with Blasius' equation.

    Parameters
    ----------
    x_array : np.ndarray
        Array of x coordinates along the airfoil.
    velocity_array : np.ndarray
        Array of external velocity along the airfoil.
    rho : float
        Density of the fluid.
    nu : float
        Kinematic viscosity of the fluid.

    Returns
    -------
    float
        Linear drag.
    """

    # Compute the length of the airfoil
    length = x_array[-1]

    # Assume the external velocity is uniform
    velocity = np.mean(velocity_array)

    # Check turbulent transition
    check_turbulent_transition(length, velocity, nu)

    # Compute the drag
    linear_drag = 2 * 0.332 * rho * \
        np.power(velocity, 3 / 2) * np.sqrt(nu) * \
        np.sqrt(length)

    if return_array:
        tau_w_array = 0.332 * rho * \
            np.power(velocity_array, 2) * \
            np.sqrt(nu / (velocity_array * x_array))
        return linear_drag, tau_w_array
    return linear_drag

def compute_polhausen_linear_drag(x_array: np.ndarray, velocity_array: np.ndarray, rho: float, nu: float, return_array: bool = False):
    """
    Compute the linear drag using the Von-Karman Polhausen's method.

    Parameters
    ----------
    x_array : np.ndarray
        Array of x coordinates along the airfoil.
    velocity_array : np.ndarray
        Array of external velocity along the airfoil.
    rho : float
        Density of the fluid.
    nu : float
        Kinematic viscosity of the fluid.

    References
    ----------
    This method is explained in details in the "Aerodynamics for Engineering Students 6th edition" page 567.

    https://univ-scholarvox-com.ezproxy.universite-paris-saclay.fr/book/88809524

    Returns
    -------
    float
        Linear drag.
    """

    # Compute the velocity gradient
    velocity_derivative_array = np.gradient(velocity_array, x_array)

    # Define the functions for the ODE
    def F1(boundary_progress):
        I = (1 / 63) * ((37 / 5) - (boundary_progress / 15) -
                        (np.power(boundary_progress, 2) / 144))
        res = I - 2 * (boundary_progress / 63) * \
            ((1 / 15) + (boundary_progress / 72))
        return res

    def F2(boundary_progress):
        res = 4 + boundary_progress / 3 - \
            (3 / 10) * boundary_progress + \
            np.power(boundary_progress, 2) / 60 + (4 / 63) *\
            ((37 / 5) - boundary_progress / 15 -
             np.power(boundary_progress, 2) / 144)
        return res

    def get_boundary_progress_from_z(x, z):
        current_velocity_array = np.interp(x, x_array, velocity_array)
        current_velocity_derivative_array = np.interp(
            x, x_array, velocity_derivative_array)
        boundary_progress = z * \
            current_velocity_derivative_array / current_velocity_array
        return boundary_progress

    def get_delta_from_z(x, z):
        current_velocity_array = np.interp(x, x_array, velocity_array)
        delta = np.sqrt(z * nu / current_velocity_array)
        return delta

    def get_z_derivative(x: np.ndarray, z: np.ndarray):
        V_s = 0  # Vs=0 for now, general case might be implemented later
        boundary_progress = get_boundary_progress_from_z(x, z)
        delta = get_delta_from_z(x, z)
        F1_val = F1(boundary_progress)
        F2_val = F2(boundary_progress)
        z_derivative = (F2_val / F1_val) + boundary_progress - \
            (2 / F1_val) * ((V_s * delta) / nu)
        # z_derivative = np.array(z_derivative)
        return z_derivative

    # Should be zero but initialized slightly above to avoid divergence
    z_0 = np.array([2e-5])
    dx = (x_array[1] - x_array[0])

    # Solve the ODE with Runge-Kutta 4
    solution = solve_ivp(
        get_z_derivative,
        t_span=[0, x_array[-1]],
        y0=z_0,
        method="RK45",
        first_step=dx,
        t_eval=x_array
    )

    # Compute boundary layer variables
    boundary_progress_array = get_boundary_progress_from_z(
        solution.t, solution.y[0])
    delta_array = get_delta_from_z(solution.t, solution.y[0])

    # Raise error if boundary layer separation
    if (boundary_progress_array < - 12).any():
        raise ValueError(
            f"Boundary layer separation, unable to compute the drag.")

    # Raise error if excess velocity
    if (boundary_progress_array > 12).any():
        raise ValueError(
            f"Excess velocity detected, unable to compute the drag. Please use a more precise model.")

    # Check turbulent transition
    check_turbulent_transition(x_array, velocity_array, nu)

    # Compute the wall shear stress
    tau_w_array = (2 + boundary_progress_array / 6) * \
        rho * nu * velocity_array / delta_array

    # Integrate to determine the drag
    linear_drag = simpson(tau_w_array, solution.t)

    if return_array:
        return linear_drag, tau_w_array
    return linear_drag

def compute_simulation_linear_drag():
    raise NotImplementedError

def compute_linear_drag(x_array: np.ndarray, velocity_array: np.ndarray, rho: float, nu: float, method: Literal["blasius", "polhausen", "simulation"], return_array: bool = False):
    if method == "blasius":
        linear_drag = compute_blasius_linear_drag(
            x_array, velocity_array, rho, nu, return_array=return_array)
    elif method == "polhausen":
        linear_drag = compute_polhausen_linear_drag(
            x_array, velocity_array, rho, nu, return_array=return_array)
    elif method == "simulation":
        linear_drag = compute_simulation_linear_drag(
            x_array, velocity_array, rho, nu, return_array=return_array)
    else:
        raise NotImplementedError(
            f"The method {method} is not implemented. Please use 'blasius', 'polhausen' or 'simulation'.")

    return linear_drag

def compute_blasius_rectangle_drag(width: float, length: float, velocity: float, rho: float, nu: float):
    """
    Compute drag produced by the upper face of a rectangle with no incidence at the given velocity.

    Parameters
    ----------
    width : float
        Width of the rectangle.
    length : float
        Length of the rectangle.
    velocity : float
        Velocity of the rectangle.
    rho : float
        Density of the fluid.
    nu : float
        Kinematic viscosity of the fluid.

    Returns
    -------
    float
        Drag of the upper face of the rectangle.
    """

    # Prepare data
    x_array = np.array([0, length])
    velocity_array = np.array([velocity, velocity])

    # Compute drag
    linear_drag = compute_blasius_linear_drag(x_array, velocity_array, rho, nu)
    drag = linear_drag * width

    return drag
