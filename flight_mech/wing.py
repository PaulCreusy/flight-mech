"""
Module to analyse wing aerodynamic properties.
"""

###########
# Imports #
###########

# Python imports #

import os
from typing import Literal

# Dependencies #

import numpy as np
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt
from copy import deepcopy

# Local imports #

from flight_mech.airfoil import Airfoil

#############
# Constants #
#############

# Define a default plane database location
default_wing_database = os.path.join(
    os.path.dirname(__file__), "wing_database")

###########
# Classes #
###########

class Wing:
    _y_array: np.ndarray
    twisting_angle_array: np.ndarray
    name: str
    chord_length_array: np.ndarray
    x_center_offset_array: np.ndarray
    base_airfoil: Airfoil

    @property
    def y_array(self) -> np.ndarray:
        """
        Array of y coordinates used for the wing definition.
        """
        return self._y_array

    @y_array.setter
    def y_array(self, value: np.ndarray):
        if self.twisting_angle_array is not None and self.chord_length_array is not None and self.x_offset_array is not None:
            if self._y_array.size == self.chord_length_array.size:
                self.re_interpolate(value, update_y_array=False)
            self._y_array = value
        else:
            self._y_array = value
        self._chord_length = np.max(value)

    def re_interpolate(self, new_y_array: np.ndarray, update_y_array: bool = True):
        pass

    def plot_2D(self):
        """
        Plot the shape of the wing in 2D.
        """
        pass

    def plot_3D(self, nb_points_airfoil=50):
        """
        Plot the shape of the wing in 3D.
        """

        # Create a display airfoil with normalized chord and less points
        display_airfoil = deepcopy(self.base_airfoil)
        display_airfoil.chord_length = 1
        display_airfoil.re_interpolate_with_cosine_distribution(
            nb_points_airfoil)

        # Create an array containing all the points
        points_array = np.zeros((nb_points_airfoil * self.y_array.size * 2, 3))
        for i in range(self.y_array.size):
            ratio = self.chord_length_array[i] / display_airfoil.chord_length
            points_array[i:i + nb_points_airfoil, 0] = display_airfoil.x_selig_array * \
                ratio + self.x_center_offset_array[i]
            points_array[i:i + nb_points_airfoil,
                         1] = display_airfoil.z_selig_array
            # TODO
