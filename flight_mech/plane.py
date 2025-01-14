"""
Module to compute the plane performances.
"""

###########
# Imports #
###########

# Python imports #

import os
from typing import Literal
import json

# Dependencies #

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Local imports #

from .atmosphere import (
    compute_air_density_from_altitude,
    compute_altitude_from_sigma,
    compute_sigma_from_altitude
)

#############
# Constants #
#############

# Define physical constants
P0 = 1e5  # Pa
g = 9.81  # m.s-2

# Define a default plane database location
default_plane_database = os.path.join(
    os.path.dirname(__file__), "plane_database")

###########
# Classes #
###########

class Plane:
    """
    Model to define a plane and compute its characteristics.
    """

    # Mass
    m_empty: float | None = None  # kg
    m_fuel: float = 0.  # kg
    m_payload: float = 0.  # kg
    P: float | None = None  # N

    # Geometry
    S: float | None = None  # m2
    b: float | None = None  # m
    c: float | None = None  # m
    wing_shape_coefficient: float = 1.
    wing_to_ground_height: float | None = None  # m

    # Thrust
    nb_engines: int = 1
    thrust_per_engine: float | None = None  # N
    engine_type: Literal["moto-propeller", "turbo-reactor"] = "turbo-reactor"
    fuel_specific_conso: float | None = None  # kg/(N.h)

    # Aero coefficients
    C_D_0: float | None = None
    k: float | None = None
    C_L_alpha: float | None = None  # rad-1
    C_L_delta: float | None = None  # rad-1
    C_m_0: float | None = None
    C_m_alpha: float | None = None  # rad-1
    C_m_delta: float | None = None  # rad-1
    a: float | None = None
    alpha_0: float | None = 0.  # rad
    alpha_stall: float = 15 * np.pi / 180  # rad
    ground_effect_coefficient: float | None = None
    C_L_max: float | None = None

    def __init__(self,
                 plane_data_name: str | None = None,
                 plane_database_folder: str = default_plane_database):
        if plane_data_name is not None:
            self.load_plane_data(plane_data_name, plane_database_folder)

    @property
    def m(self) -> float:
        """
        Total mass of the plane.
        """

        return self.m_empty + self.m_fuel + self.m_payload

    @property
    def extension(self) -> float:
        """
        Extension coefficient of the wings.
        """

        return pow(self.b, 2) / self.S

    @property
    def f_max(self) -> float:
        """
        Max gliding ratio.
        """

        f_max = 1 / (2 * np.sqrt(self.k * self.C_D_0))

        return f_max

    @property
    def C_L_f_max(self):
        """
        Lift coefficient at max gliding ratio.
        """

        C_L_star = np.sqrt(self.C_D_0 / self.k)

        return C_L_star

    @property
    def alpha_f_max(self):
        """
        Angle of incidence at max gliding ratio.
        """

        C_L_star = self.C_L_f_max
        alpha_f_max = C_L_star / self.a + self.alpha_0

        return alpha_f_max

    @property
    def fuel_specific_conso_SI(self):
        """
        Fuel specific consumption in International unit system.
        """

        return self.fuel_specific_conso / 3600

    def update_k(self, force: bool = False):
        """
        Update the value of k.

        Parameters
        ----------
        force : bool, optional
            Force the computation of k even if it already exists, by default False
        """

        if self.k is None or force:
            self.k = 1 / (np.pi * self.wing_shape_coefficient * self.extension)

    def update_P(self, force: bool = False):
        """
        Update the value of the weight.

        Parameters
        ----------
        force : bool, optional
            Force the computation of the weight even if it already exists, by default False
        """

        if self.P is None or force:
            self.P = self.m * g

    def update_ground_effect_coefficient(self, force: bool = False):
        if (self.ground_effect_coefficient is None or force) and self.wing_to_ground_height is not None:
            temp = np.power(16 * self.wing_to_ground_height / self.b, 2)
            self.ground_effect_coefficient = temp / (1 + temp)

    def update_C_L_max(self, force: bool = False):
        if (self.C_L_max is None or force) and\
                self.alpha_stall is not None and self.a is not None:
            self.C_L_max = self.C_L(self.alpha_stall)

    def update_variables(self, force: bool = False):
        self.update_k(force)
        self.update_P(force)
        self.update_ground_effect_coefficient(force)
        self.update_C_L_max(force)

    def C_L(self, alpha: float):
        C_L = self.a * (alpha - self.alpha_0)
        return C_L

    def C_D(self, alpha: float = None, C_L: float = None):
        self.update_k()
        if C_L is None:
            C_L = self.C_L(alpha)
        C_D = self.C_D_0 + self.k * pow(C_L, 2)
        return C_D

    def f(self, alpha: float):
        C_L = self.C_L(alpha)
        C_D = self.C_D(alpha)
        f = C_L / C_D
        return f

    def v(self, alpha: float, z: float = 0.):
        rho = compute_air_density_from_altitude(z)
        v = np.sqrt(self.P / (.5 * rho * self.S * self.C_L(alpha)))
        return v

    def n(self, v: float, z: float, alpha: float):
        lift = self.compute_lift(v, z, alpha)
        n = lift / self.P
        return n

    def compute_drag(self,
                     v: float,
                     z: float,
                     alpha: float | None = None,
                     C_L: float | None = None):
        drag = .5 * \
            compute_air_density_from_altitude(
                z) * self.S * pow(v, 2) * self.C_D(alpha, C_L)
        return drag

    def compute_lift(self,
                     v: float,
                     z: float,
                     alpha: float | None = None,
                     C_L: float | None = None):
        if C_L is None:
            C_L = self.C_L(alpha)
        lift = .5 * \
            compute_air_density_from_altitude(z) * self.S * pow(v, 2) * C_L
        return lift

    def compute_thrust(self, z: float):
        if self.engine_type == "turbo-reactor":
            sigma = compute_sigma_from_altitude(z)
            return self.thrust_per_engine * self.nb_engines * sigma
        else:
            raise NotImplementedError

    def compute_normalized_thrust(self, z: float):
        t = self.compute_thrust(z) * self.f_max / self.P
        return t

    def compute_stall_speed(self,
                            z: float = 0.,
                            alpha: float | None = None,
                            C_L_max: float | None = None):
        if alpha is None:
            alpha = self.alpha_stall
        if C_L_max is None:
            C_L_max = self.C_L(alpha)
        rho = compute_air_density_from_altitude(z)
        stall_speed = np.sqrt((2 * self.P) / (rho * self.S * C_L_max))
        return stall_speed

    def compute_gliding_speed(self, alpha: float, z: float):
        rho = compute_air_density_from_altitude(z)
        v = np.sqrt((2 * self.P) / (rho * self.S * self.C_L(alpha)))
        return v

    def compute_min_descent_gliding_slope(self):
        gamma_min = -1 / self.f_max
        return gamma_min

    def compute_v_at_gliding_v_z_min(self, z: float):
        rho = compute_air_density_from_altitude(z)
        v_at_v_z_min = np.sqrt((2 * self.P) / (rho * self.S)) * \
            pow(self.k / (3 * self.C_D_0), 1 / 4)
        return v_at_v_z_min

    def compute_gliding_v_z_min(self, z: float):
        rho = compute_air_density_from_altitude(z)
        v_z_min = 4 * np.sqrt((2 * self.P) / (rho * self.S)) * \
            pow(pow(self.k, 3) * self.C_D_0 / 27, 1 / 4)
        return v_z_min

    def compute_max_gliding_time(self, z: float):
        v_z_min = self.compute_gliding_v_z_min(z)
        return z / v_z_min

    def compute_velocity_interval_for_fixed_thrust(self, z: float):
        t = self.compute_normalized_thrust(z)
        v_ref = self.compute_reference_speed(z)
        v_max = float(np.sqrt(t + np.sqrt(pow(t, 2) - 1)) * v_ref)
        v_min = float(np.sqrt(t - np.sqrt(pow(t, 2) - 1)) * v_ref)
        return v_min, v_max

    def compute_reference_speed(self, z: float):
        """
        Corresponds to the gliding velocity at f_max.

        Parameters
        ----------
        z : float
            Altitude of the plane.

        Returns
        -------
        float
            Reference velocity.
        """

        rho = compute_air_density_from_altitude(z)
        v_ref = np.sqrt((2 * self.P) / (rho * self.S)) * \
            pow(self.k / self.C_D_0, 1 / 4)

        return v_ref

    def compute_max_gliding_range(self, z: float):
        R_max = self.f_max * z
        return R_max

    def compute_thrust_needed(self, alpha: float, z: float):
        rho = compute_air_density_from_altitude(z)
        friction_drag = .5 * rho * self.S * self.v(alpha, z) * self.C_D_0
        induced_drag = (2 * self.k * self.P) / \
            (rho * self.S * pow(self.v(alpha, z), 2))
        thrust_needed = friction_drag + induced_drag
        return thrust_needed

    def compute_min_thrust_needed(self):
        min_thrust_needed = 2 * self.P * np.sqrt(self.k * self.C_D_0)
        return min_thrust_needed

    def compute_speed_for_min_thrust_needed(self, z: float):
        return self.compute_reference_speed(z)

    def compute_speed_for_min_power_needed(self, z: float):
        reference_speed = self.compute_reference_speed(z)
        v_for_w_min = (1 / pow(3, 1 / 4)) * reference_speed
        return v_for_w_min

    def compute_max_altitude(self):
        thrust_needed_at_f_max = self.P / self.f_max
        sigma = thrust_needed_at_f_max / self.compute_thrust(0)
        z_max = compute_altitude_from_sigma(sigma)
        return z_max

    def compute_speed_for_max_ascension_speed(self, z: float):
        t = self.compute_normalized_thrust(z)
        u_m = np.sqrt((t + np.sqrt(np.power(t, 2) + 3)) / 3)
        v_ref = self.compute_reference_speed(z)
        return v_ref * u_m

    def compute_max_ascension_speed(self, z: float):
        t = self.compute_normalized_thrust(z)
        u_m = np.sqrt((t + np.sqrt(np.power(t, 2) + 3)) / 3)
        v_ref = self.compute_reference_speed(z)
        v_z_max = (1 / (2 * self.f_max)) * (2 * t * u_m -
                                            (np.power(u_m, 3) + 1 / u_m)) * v_ref
        return v_z_max

    def compute_ascension_slope(self, alpha: float, z: float):
        T = self.compute_thrust(z)
        gamma = np.arcsin(T / self.P - 1 / self.f(alpha))
        return gamma

    def compute_max_ascension_slope(self, z: float):
        t = self.compute_normalized_thrust(z)
        gamma_max = np.arcsin((t - 1) / self.f_max)
        return gamma_max

    def compute_load_factor_from_roll_angle(self, phi: float):
        n_z = 1 / np.cos(phi)
        return n_z

    def compute_max_range_at_fixed_altitude(self, z: float):
        rho = compute_air_density_from_altitude(z)
        optimal_C_L = self.C_L_f_max / np.sqrt(3)
        plane_range = (2 / (self.fuel_specific_conso_SI * g))\
            * (np.sqrt(optimal_C_L) / self.C_D(C_L=optimal_C_L)) * \
            np.sqrt(2 / (rho * self.S)) * \
            (np.sqrt(self.P) - np.sqrt((self.m_empty + self.m_payload) * g))
        return plane_range

    def compute_range_at_fixed_speed(self,
                                     v: float,
                                     alpha: float | None = None,
                                     f: float | None = None):
        if f is None:
            f = self.f(alpha)
        plane_range = (v / (self.fuel_specific_conso_SI * g)) * \
            f * np.log(self.P / ((self.m_empty + self.m_payload) * g))
        return plane_range

    def compute_endurance(self,
                          alpha: float | None = None,
                          f: float | None = None):
        if f is None:
            f = self.f(alpha)
        endurance = (1 / (self.fuel_specific_conso_SI * g)) * f * \
            np.log(self.P / ((self.m_empty + self.m_payload) * g))
        return endurance

    def compute_take_off_distance_no_friction(self, z: float):
        T = self.compute_thrust(z)
        rho = compute_air_density_from_altitude(z)
        d_take_off = (self.P / T) * (1.44 * self.P / self.S) / \
            (rho * g * self.C_L_max)
        return d_take_off

    def compute_ground_effect(self,
                              alpha: float | None = None,
                              C_L: float | None = None):
        if C_L is None:
            C_L = self.C_L(alpha)
        ground_effect = self.ground_effect_coefficient * \
            np.power(C_L, 2) / (self.wing_shape_coefficient *
                                np.pi * self.extension)
        return ground_effect

    def compute_drag_with_ground_effect(self,
                                        v: float,
                                        z: float,
                                        alpha: float | None = None,
                                        C_L: float | None = None):
        ground_effect = self.compute_ground_effect(alpha=alpha, C_L=C_L)
        rho = compute_air_density_from_altitude(z)
        drag_with_ground_effect = .5 * rho * self.S * \
            np.power(v, 2) * (self.C_D_0 + ground_effect)
        return drag_with_ground_effect

    def compute_take_off_speed(self, z: float):
        rho = compute_air_density_from_altitude(z)
        take_off_speed = 1.2 * \
            np.sqrt(2 * self.P / (rho * self.S * self.C_L_max))
        return take_off_speed

    def compute_landing_speed(self, z: float):
        rho = compute_air_density_from_altitude(z)
        take_off_speed = 1.3 * \
            np.sqrt(2 * self.P / (rho * self.S * self.C_L_max))
        return take_off_speed

    def compute_take_off_distance_with_friction(self,
                                                z: float,
                                                mu: float,
                                                C_L_max: float | None = None):
        if C_L_max is None:
            C_L_max = self.C_L_max
        take_off_speed = self.compute_take_off_speed(z)
        T = self.compute_thrust(z)
        D = self.compute_drag_with_ground_effect(
            take_off_speed * 0.7, z, C_L=C_L_max)
        L = self.compute_lift(take_off_speed * 0.7, z, C_L=C_L_max)
        rho = compute_air_density_from_altitude(z)
        d_take_off = (1.44 * (self.P / self.S)) / (rho * g *
                                                   C_L_max) * (self.P / (T - (D + mu * (self.P - L))))
        return d_take_off

    def compute_landing_distance(self,
                                 z: float,
                                 mu: float,
                                 reverse_thrust: float = 0.,
                                 C_L: float | None = None,
                                 C_L_max: float | None = None):
        if C_L_max is None:
            C_L_max = self.C_L_max
        if C_L is None:
            C_L = C_L_max

        landing_speed = self.compute_landing_speed(z)
        D = self.compute_drag(landing_speed * 0.7, z, C_L=C_L)
        L = self.compute_lift(landing_speed * 0.7, z, C_L=C_L)
        rho = compute_air_density_from_altitude(z)
        d_landing = (1.69 * (self.P / self.S)) / (rho * g * C_L_max) * \
            (self.P / (reverse_thrust + (D + mu * (self.P - L))))

        return d_landing

    def compute_alpha_and_delta_at_flight_point(self,
                                                z: float = 0,
                                                v: float | None = None):
        rho = compute_air_density_from_altitude(z)
        if v is None:
            v = self.compute_reference_speed(z)
        A = np.array([
            [self.C_L_alpha, self.C_L_delta],
            [self.C_m_alpha, self.C_m_delta]
        ])
        B = np.array([
            [self.P / (.5 * rho * self.S * np.power(v, 2))],
            [-self.C_m_0]
        ])
        sol = np.linalg.solve(A, B)
        alpha = self.alpha_0 + sol[0]
        delta = sol[1]
        return alpha, delta

    def plot_polar_graph(self, nb_points: int = 100):
        alpha_array = np.linspace(self.alpha_0, self.alpha_stall, nb_points)
        C_L_array = self.C_L(alpha_array)
        C_D_array = self.C_D(alpha_array)
        plt.plot(C_D_array, C_L_array)
        plt.xlabel("C_D")
        plt.ylabel("C_L")
        plt.title("Polar graph of the plane")
        plt.show()

    def plot_gliding_TV_graph(self, z: float | list | tuple = 0., nb_points: int = 100):
        colors_list = list(mcolors.TABLEAU_COLORS.values())
        if not isinstance(z, (list, tuple)):
            z = [z]
        for i in range(len(z)):
            alpha_array = np.linspace(
                self.alpha_0, self.alpha_stall, nb_points)
            T_array = self.P / self.f(alpha_array)
            V_array = self.v(alpha_array, z[i])
            plt.plot(V_array, T_array, label=f"z={z[i]}", color=colors_list[i])
            if self.engine_type == "turbo-reactor" and self.thrust_per_engine is not None:
                thrust = self.compute_thrust(z[i])
                plt.plot([0, np.max(V_array[V_array != np.inf])], [
                         thrust, thrust], "--", label=f"Thrust at z={z[i]}", color=colors_list[i])
        plt.legend()
        plt.xlabel("v")
        plt.ylabel("T")
        plt.title("TV graph of the plane")

        plt.show()

    def plot_gliding_WV_graph(self, z: float | list | tuple = 0., nb_points: int = 100):
        if not isinstance(z, (list, tuple)):
            z = [z]
        for i in range(len(z)):
            alpha_array = np.linspace(
                self.alpha_0, self.alpha_stall, nb_points)
            T_array = self.P / self.f(alpha_array)
            V_array = self.v(alpha_array, z[i])
            W_array = T_array * V_array
            plt.plot(V_array, W_array, label=f"z={z[i]}")
        plt.xlabel("v")
        plt.ylabel("W")
        plt.title("WV graph of the plane")
        plt.show()

    def load_plane_data(self, plane_data_name: str, plane_data_folder: str = default_plane_database):
        """
        Load the data from a plane stored in the given database folder.

        Parameters
        ----------
        plane_data_name : str
            Name of the plane.
        plane_data_folder : str, optional
            Path to the database folder, by default default_plane_database
        """

        # Load the json file
        file_path = os.path.join(plane_data_folder, plane_data_name + ".json")
        with open(file_path, "r") as file:
            plane_data_dict = json.load(file)

        # List the variables allowed in the class
        plane_allowed_variables_list = dir(self)

        # Store the values in the class
        for key in plane_data_dict:
            if key in plane_allowed_variables_list:
                setattr(self, key, plane_data_dict[key])

        # Update the variables
        self.update_variables()
