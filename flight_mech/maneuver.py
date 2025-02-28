"""
Module to define plane maneuvers.
"""

###########
# Imports #
###########

# Python imports #

from typing import Literal, Callable
from abc import ABC, abstractmethod

# Dependencies #

import numpy as np

# Local imports #

from flight_mech.plane import (
    Plane
)
from flight_mech._common import plot_graph

###########
# Classes #
###########

class Maneuver(ABC):
    """
    Base class to create maneuvers.
    """

    time_list: list[float]
    evolution_variable_names: list[str]
    evolution_variable_units: dict[str, str]

    def plot_graph(self, variable: str, **kwargs):
        """
        Plot the graph of evolution of the given variable during the maneuver.

        Parameters
        ----------
        variable : str
            Name of the variable to plot.

        Note
        ----
        For more details on the optional arguments, please check flight_mech._common.plot_graph.
        """

        # Plot
        plot_graph(
            x_array=self.time_list,
            y_array=self.__getattribute__(f"{variable}_list"),
            title=f"{variable.capitalize()} graph",
            x_label=f"Time [{self.evolution_variable_units['time']}]",
            y_label=f"{variable.capitalize()} [{self.evolution_variable_units[variable]}]",
            **kwargs
        )

    @abstractmethod
    def compute_evolution(self) -> None:
        """
        Compute the evolution of variables during the maneuver.
        """
        pass

    @abstractmethod
    def _initialize_evolution_variables(self):
        """
        Initialize the evolution variables.
        """
        pass

class TakeOffManeuver(Maneuver):
    """
    Class to compute the evolution of variables during a take off sequence.
    """

    thrust_evolution_function: Callable[[float], float]
    ground_rotation_speed: float  # rad.s-1
    air_rotation_speed: float = 1.0  # rad.s-1
    rotation_sequence_trigger_speed: float  # m.s-1
    ground_friction_coefficient: float
    max_ground_angle_of_incidence: float  # rad
    end_take_off_altitude: float  # m
    plane_model: Plane
    dt: float = 0.1
    evolution_variable_names = [
        "phase",
        "altitude",
        "time",
        "thrust",
        "incidence",
        "pitch",
        "lift_coefficient",
        "drag_coefficient",
        "x_coord",
        "velocity",
        "acceleration",
        "ground_reaction_force",
        "ground_friction_force",
        "lift",
        "drag",
        "rho"
    ]
    _nb_max_iterations: int = 1000
    rotation_start_time: float | None = None
    flight_start_time: float | None = None
    evolution_variable_units = {
        "phase": "",
        "altitude": "m",
        "time": "s",
        "thrust": "N",
        "incidence": "rad",
        "pitch": "rad",
        "lift_coefficient": "",
        "drag_coefficient": "",
        "x_coord": "m",
        "velocity": "m.s-1",
        "acceleration": "m.s-2",
        "ground_reaction_force": "N",
        "ground_friction_force": "N",
        "lift": "N",
        "drag": "N",
        "rho": "kg.m-3"
    }

    def __init__(self,
                 plane_model: Plane,
                 rotation_speed: float,
                 rotation_sequence_trigger_speed: float,
                 thrust_evolution_function: Callable[[
                     float], float] = lambda x: 1.,
                 ground_friction_coefficient: float = 0.03,
                 max_ground_angle_of_incidence: float = 12 * np.pi / 180,
                 ground_altitude: float = 0.,
                 end_take_off_altitude: float | None = None):

        self.plane_model = plane_model
        self.ground_rotation_speed = rotation_speed
        self.thrust_evolution_function = thrust_evolution_function
        self.rotation_sequence_trigger_speed = rotation_sequence_trigger_speed
        self.ground_friction_coefficient = ground_friction_coefficient
        self.max_ground_angle_of_incidence = max_ground_angle_of_incidence
        self.ground_altitude = ground_altitude
        if end_take_off_altitude is None:
            self.end_take_off_altitude = ground_altitude + 30

    def _initialize_evolution_variables(self):
        self.phase_list = ["initial_acceleration"]
        self.time_list = [0.]
        self.thrust_list = [
            self.thrust_evolution_function(0)]
        self.incidence_list = [0.]
        self.pitch_list = [0.]
        self.altitude_list = [self.ground_altitude]
        self.lift_coefficient_list = [
            self.plane_model.C_L(0)]
        self.drag_coefficient_list = [self.plane_model.C_D_0]
        self.x_coord_list = [0.]
        self.velocity_list = [0.]
        self.acceleration_list = [0.]
        self.pitch_derivative_list = [0.]
        self.ground_reaction_force_list = [
            self.plane_model.P]
        self.ground_friction_force_list = [
            self.plane_model.P * self.ground_friction_coefficient]
        self.drag_list = [self.plane_model.compute_drag(
            0, z=self.ground_altitude, alpha=0)]
        self.lift_list = [self.plane_model.compute_lift(
            0, z=self.ground_altitude, alpha=0)]
        self.rho_list = [self.plane_model.atmosphere_model.compute_density_from_altitude(
            self.ground_altitude)]

        # Reset timers
        self.rotation_start_time = None
        self.flight_start_time = None

    def compute_evolution(self):

        # Initialise
        self._initialize_evolution_variables()

        # Iterate until lift off
        while self.altitude_list[-1] < self.end_take_off_altitude and len(self.phase_list) < self._nb_max_iterations:
            # Compute time
            current_time = self.time_list[-1] + self.dt

            # Determine the phase
            if self.velocity_list[-1] > self.rotation_sequence_trigger_speed and self.ground_reaction_force_list[-1] <= 0:
                current_phase = "flight"
                if self.phase_list[-1] == "rotation":
                    self.flight_start_time = current_time
            elif self.velocity_list[-1] > self.rotation_sequence_trigger_speed:
                current_phase = "rotation"
                if self.phase_list[-1] == "initial_acceleration":
                    self.rotation_start_time = current_time
            else:
                current_phase = "initial_acceleration"

            # Compute thrust
            current_thrust = self.plane_model.thrust_per_engine * \
                self.plane_model.nb_engines * \
                self.thrust_evolution_function(current_time)

            # Compute angle of incidence
            current_incidence = min(self.incidence_list[-1] + self.ground_rotation_speed * self.dt * (
                current_phase == "rotation") + self.air_rotation_speed * self.dt * (current_phase == "flight"), self.plane_model.alpha_stall)

            # Compute pitch
            if current_phase == "flight":
                current_pitch = self.pitch_list[-1] + \
                    self.pitch_derivative_list[-1] * self.dt
            else:
                current_pitch = 0

            # Compute lift coefficient
            current_lift_coefficient = self.plane_model.C_L(current_incidence)

            # Compute drag coefficient
            current_drag_coefficient = self.plane_model.C_D(current_incidence)

            # Compute rho
            current_rho = self.plane_model.atmosphere_model.compute_density_from_altitude(
                self.altitude_list[-1])

            # Compute drag
            current_drag = self.plane_model.compute_drag(
                v=self.velocity_list[-1],
                z=self.altitude_list[-1],
                alpha=self.incidence_list[-1]
            )

            # Compute lift
            current_lift = self.plane_model.compute_lift(
                v=self.velocity_list[-1],
                z=self.altitude_list[-1],
                alpha=self.incidence_list[-1]
            )

            # Compute pitch derivative
            vertical_force = -(-self.plane_model.P * np.cos(current_pitch) +
                               current_lift + current_thrust * np.sin(current_incidence))
            if current_phase == "flight":
                current_pitch_derivative = -vertical_force / \
                    (self.plane_model.m * self.velocity_list[-1])
            else:
                current_pitch_derivative = 0

            # Compute ground reaction force
            if current_phase == "flight":
                current_ground_reaction_force = 0
            else:
                current_ground_reaction_force = max(vertical_force, 0)

            # Compute ground friction force
            current_ground_friction_force = current_ground_reaction_force * \
                self.ground_friction_coefficient

            # Compute acceleration
            if current_phase == "flight":
                current_acceleration = (current_thrust * np.cos(current_incidence) -
                                        current_drag - self.plane_model.P * np.sin(current_pitch)) / self.plane_model.m
            else:
                current_acceleration = (current_thrust * np.cos(current_incidence) -
                                        current_drag - current_ground_friction_force) / self.plane_model.m

            # Compute velocity
            current_velocity = self.velocity_list[-1] + \
                current_acceleration * self.dt

            # Compute distance
            current_x_coord = self.x_coord_list[-1] + \
                current_velocity * self.dt * np.cos(current_pitch)

            # Compute altitude
            current_altitude = self.altitude_list[-1] + \
                current_velocity * self.dt * np.sin(current_pitch)

            # Add the variables in the lists
            self.phase_list.append(current_phase)
            self.time_list.append(current_time)
            self.thrust_list.append(current_thrust)
            self.incidence_list.append(current_incidence)
            self.pitch_list.append(current_pitch)
            self.lift_coefficient_list.append(current_lift_coefficient)
            self.drag_coefficient_list.append(current_drag_coefficient)
            self.acceleration_list.append(current_acceleration)
            self.velocity_list.append(current_velocity)
            self.x_coord_list.append(current_x_coord)
            self.rho_list.append(current_rho)
            self.pitch_derivative_list.append(current_pitch_derivative)
            self.drag_list.append(current_drag)
            self.lift_list.append(current_lift)
            self.altitude_list.append(current_altitude)
            self.ground_reaction_force_list.append(
                current_ground_reaction_force)
            self.ground_friction_force_list.append(
                current_ground_friction_force)
