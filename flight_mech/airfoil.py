"""
Module to analyse airfoil aerodynamic properties.
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
import pandas as pd

###########
# Classes #
###########

class Airfoil:

    extrados_y_array: np.ndarray
    intrados_y_array: np.ndarray
    x_array: np.ndarray

    def __init__(self):
        pass

    @property
    def camber_y_array(self) -> np.ndarray:
        camber_y_array = (self.extrados_y_array + self.intrados_y_array) / 2
        return camber_y_array

    @property
    def max_thickness(self):
        pass

    @property
    def max_camber(self):
        pass

    def load_selig_file(self, file_path: str, skiprows: int = 0):
        file_content = np.loadtxt(file_path, skiprows=skiprows)
        x_selig_array = file_content[:, 0]
        y_selig_array = file_content[:, 1]
        self.import_xy_selig_arrays(x_selig_array, y_selig_array)

    def import_xy_selig_arrays(self, x_selig_array: np.ndarray, y_selig_array: np.ndarray):
        # Compute the x derivative to split parts of the airfoil
        x_diff = np.zeros(x_selig_array.shape)
        x_diff[:-1] = x_selig_array[1:] - x_selig_array[:-1]
        x_diff[-1] = x_diff[-2]

        # Split parts
        part_1_x = x_selig_array[x_diff * x_diff[0] > 0]
        part_1_y = y_selig_array[x_diff * x_diff[0] > 0]
        part_2_x = x_selig_array[x_diff * x_diff[0] < 0]
        part_2_y = y_selig_array[x_diff * x_diff[0] < 0]

        # Extract all x locations
        self.x_array = np.unique(x_selig_array)

        # Reorder arrays for interpolation
        part_1_order = np.argsort(part_1_x)
        part_2_order = np.argsort(part_2_x)

        # Interpolate both parts on all x locations
        part_1_y_interpolated = np.interp(
            self.x_array, part_1_x[part_1_order], part_1_y[part_1_order])
        part_2_y_interpolated = np.interp(
            self.x_array, part_2_x[part_2_order], part_2_y[part_2_order])

        # Assign depending on which part is extrados or intrados
        if np.mean(part_1_y_interpolated) > np.mean(part_2_y_interpolated):
            self.extrados_y_array = part_1_y_interpolated
            self.intrados_y_array = part_2_y_interpolated
        else:
            self.extrados_y_array = part_2_y_interpolated
            self.intrados_y_array = part_1_y_interpolated

    def plot(self):
        x_plot = np.concatenate(
            (self.x_array, self.x_array[::-1], [self.x_array[0]]), axis=0)
        y_plot = np.concatenate(
            (self.extrados_y_array, self.intrados_y_array[::-1], [self.extrados_y_array[0]]), axis=0)
        plt.plot(x_plot, y_plot)
        plt.axis("equal")
        plt.show()
