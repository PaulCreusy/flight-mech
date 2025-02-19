"""
Module to define motor models.
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
from scipy.optimize import minimize

# Local imports #

from flight_mech._common import plot_graph

#############
# Constants #
#############

# Define a default plane database location
default_motor_database = os.path.join(
    os.path.dirname(__file__), "motor_database")

# Define the units of the variables to indicate them on the plot labels
VARIABLE_TO_UNIT = {
    "efficiency": ".",
    "rotation_speed": "rad.s-1",
    "power": "W",
    "torque": "N.m"
}

###########
# Classes #
###########

class ElectricalMotor:
    """
    Class to define an electrical motor controlled by tension.

    Parameters
    ----------
    I0 : float
        No-load current in A.
    internal_resistance : float
        Internal resistance of the motor in Ohm.
    external_resistance : float
        External resistance of the batteries in Ohm.
    K : float
        Coefficient between tension and rotation speed in s-1.V-1.
    """

    I_0: float  # A
    internal_resistance: float  # Ohm
    external_resistance: float  # Ohm
    Kv: float  # rad.s-1.V-1
    I_max: float  # A
    U_max: float  # V
    mass: float | None  # kg

    # Internal variables
    _U: float | None = None  # V
    _I: float | None = None  # A
    _rotation_speed: float | None = None  # rad.s-1
    _power: float | None = None  # W
    _torque: float | None = None  # N.m

    def __init__(self, motor_data_name: str | None = None, motor_database_folder: str = default_motor_database, motor_parameters_dict: dict | None = None):
        if motor_data_name is not None and motor_parameters_dict is not None:
            raise ValueError(
                "You have provided both a data file to load and a parameters dictionary. Please provide only one of them to create a motor instance.")
        elif motor_data_name is not None:
            self.load_motor_data(motor_data_name, motor_database_folder)
        elif motor_parameters_dict is not None:
            self.set_motor_parameters(**motor_parameters_dict)

    def load_motor_data(self, motor_data_name: str, motor_data_folder: str = default_motor_database):
        """
        Load the data from a motor stored in the given database folder.

        Parameters
        ----------
        motor_data_name : str
            Name of the motor.
        motor_data_folder : str, optional
            Path to the database folder, by default default_motor_database
        """

        # Load the json file
        file_path = os.path.join(motor_data_folder, motor_data_name + ".json")
        with open(file_path, "r") as file:
            motor_parameters_dict = json.load(file)

        # Load the parameters in the instance
        self.set_motor_parameters(**motor_parameters_dict)

    def set_motor_parameters(self, I_0: float, I_max: float, internal_resistance: float, external_resistance: float, Kv: float, U_max: float, mass: float | None = None):
        """
        Set the parameters of the motor.

        Parameters
        ----------
        I0 : float
            No-load current in A.
        internal_resistance : float
            Internal resistance of the motor in Ohm.
        external_resistance : float
            External resistance of the batteries in Ohm.
        Kv : float
            Coefficient between tension and rotation speed in s-1.V-1.
        U_max : float
            Max tension of the motor in V.
        mass : float | None, optional
            Mass of the motor in kg.
        """

        self.I_0 = I_0
        self.I_max = I_max
        self.internal_resistance = internal_resistance
        self.external_resistance = external_resistance
        self.Kv = Kv
        self.U_max = U_max
        self.mass = mass

    @property
    def U(self):
        if self._U is not None:
            return self._U
        else:
            return self.U_max

    @U.setter
    def U(self, value):
        if value > self.U_max:
            raise ValueError(
                f"The given tension ({value}) is above the motor maximum tension ({self.U_max}). Please provide a lower value or change the value of U_max.")
        self._U = value

    @property
    def _I_minimize_start(self):
        return (self.I_0 + self.I_max) / 2

    @property
    def I(self):
        if self._I is not None:
            return self._I
        elif self._rotation_speed is not None:
            return self._compute_I_from_rotation_speed()
        elif self._power is not None:
            return self._compute_I_from_power()
        elif self._torque is not None:
            return self._compute_I_from_torque()
        else:
            raise ValueError(
                "No value to compute the current, please specify a desired current, rotation speed, power or torque.")

    @I.setter
    def I(self, value: float):
        if value > self.I_max:
            raise ValueError(
                f"The given current ({value} A) exceeds the max current of the motor ({self.I_max} A).")
        elif value < self.I_0:
            raise ValueError(
                f"The given current ({value} A) is under the no-load current of the motor ({self.I_0} A).")
        self._I = value
        self._rotation_speed = None
        self._power = None
        self._torque = None

    @property
    def rotation_speed(self):
        if self._rotation_speed is not None:
            return self._rotation_speed
        else:
            return self._compute_rotation_speed_from_I(self.I)

    @rotation_speed.setter
    def rotation_speed(self, value: float):
        self._rotation_speed = value
        self._I = None
        self._power = None
        self._torque = None

    @property
    def power(self):
        if self._power is not None:
            return self._power
        else:
            return self._compute_power_from_I(self.I)

    @power.setter
    def power(self, value: float):
        self._power = value
        self._I = None
        self._rotation_speed = None
        self._torque = None

    @property
    def torque(self):
        if self._torque is not None:
            return self._torque
        else:
            return self._compute_torque_from_I(self.I)

    @torque.setter
    def torque(self, value: float):
        self._torque = value
        self._I = None
        self._rotation_speed = None
        self._power = None

    @property
    def efficiency(self):
        efficiency = self._compute_efficiency_from_I(self.I)
        return efficiency

    @property
    def electromotive_force(self):
        electromotive_force = self._compute_electromotive_force_from_I(self.I)
        return electromotive_force

    @property
    def max_efficiency(self):
        max_efficiency = self._compute_efficiency_from_I(
            self.compute_I_at_max_efficiency())
        return max_efficiency

    def _compute_efficiency_from_I(self, I: float):
        """
        Compute the efficiency at the given current.

        Parameters
        ----------
        I : float
            Current in the motor.

        Returns
        -------
        float
            Efficiency of the motor.
        """

        efficiency = 1 - (self.I_0 * (1 / I - self.external_resistance / self.U) + (
            (self.internal_resistance + self.external_resistance) / self.U) * I)

        return efficiency

    def _compute_electromotive_force_from_I(self, I: float):
        """
        Compute the electromotive force at the given current.

        Parameters
        ----------
        I : float
            Current in the motor.

        Returns
        -------
        float
            Electromotive force of the motor.
        """

        electromotive_force = self.U - self.internal_resistance * \
            (I - self.I_0) - self.external_resistance * I

        return electromotive_force

    def _compute_rotation_speed_from_I(self, I: float):
        """
        Compute the rotation speed at the given current.

        Parameters
        ----------
        I : float
            Current in the motor.

        Returns
        -------
        float
            Rotation speed of the motor in rad.s-1.
        """

        electromotive_force = self._compute_electromotive_force_from_I(I)
        rotation_speed = self.Kv * electromotive_force

        return rotation_speed

    def _compute_power_from_I(self, I: float):
        """
        Compute the power at the given current.

        Parameters
        ----------
        I : float
            Current in the motor.

        Returns
        -------
        float
            Power of the motor.
        """

        power = self.U * I * self._compute_efficiency_from_I(I)

        return power

    def _compute_torque_from_I(self, I: float):
        """
        Compute the torque at the given current.

        Parameters
        ----------
        I : float
            Current in the motor.

        Returns
        -------
        float
            Torque of the motor.
        """

        torque = self._compute_power_from_I(I)\
            / self._compute_rotation_speed_from_I(I)

        return torque

    def _compute_I_from_rotation_speed(self):
        """
        Compute I using the rotation speed value.

        Returns
        -------
        float
            Current in A.
        """

        I = (1 / (self.internal_resistance + self.external_resistance)) * \
            (self.U + self.internal_resistance *
             self.I_0 - self.rotation_speed / self.Kv)

        return I

    def _compute_I_from_power(self):
        """
        Compute I using the power value.

        Returns
        -------
        float
            Current in A.
        """

        a = -(self.internal_resistance + self.external_resistance) / self.U
        b = 1 + self.external_resistance * self.I_0 / self.U
        c = - self.I_0 - self.power / self.U
        delta = np.power(b, 2) - 4 * a * c
        I = (-b + np.sqrt(delta)) / (2 * a)

        return I

    def _compute_I_from_torque(self):
        """
        Compute I using the torque value.

        Returns
        -------
        float
            Current in A.
        """

        a = -(self.internal_resistance + self.external_resistance) / self.U
        b = 1 + self.external_resistance * self.I_0 / self.U + \
            (self.external_resistance + self.internal_resistance) * \
            self.torque * self.Kv / self.U
        c = - self.I_0 - (self.U + self.internal_resistance * self.I_0) * \
            self.torque * self.Kv / self.U
        delta = np.power(b, 2) - 4 * a * c
        I = (-b + np.sqrt(delta)) / (2 * a)

        return I

    def plot_graph(self, variable: Literal["efficiency", "rotation_speed", "power", "torque"], nb_points=100, **kwargs):
        """
        Plot the graph of evolution of the given variable with the current in the motor.

        Parameters
        ----------
        variable : Literal["efficiency", "rotation_speed", "power", "torque"]
            Name of the variable to plot.
        """

        # Allocate list for the x and y values
        values_array = np.zeros(nb_points)
        current_array = np.linspace(self.I_0, self.I_max, nb_points)

        # Extract the values
        for i, value in enumerate(current_array):
            values_array[i] = self.__getattribute__(
                f"_compute_{variable}_from_I")(value)

        # Plot
        plot_graph(
            x_array=current_array,
            y_array=values_array,
            title=f"{variable.capitalize()} graph",
            x_label="Current [A]",
            y_label=f"{variable.capitalize()} [{VARIABLE_TO_UNIT[variable]}]",
            **kwargs
        )

    def compute_I_at_max_efficiency(self):
        """
        Compute the current that gives the maximum efficiency.

        Returns
        -------
        float
            Current of maximum efficiency
        """

        I = np.sqrt(self.I_0 * self.U /
                    (self.internal_resistance + self.external_resistance))

        return I

    def set_at_max_efficiency(self):
        """
        Set the motor in the maximum efficiency conditions.
        """
        self.I = self.compute_I_at_max_efficiency()
