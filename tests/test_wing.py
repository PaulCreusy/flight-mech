"""
Test file for the wing module.
"""

###########
# Imports #
###########

# Python imports #

import os
import sys
sys.path.append(".")

# Dependencies #

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Local imports #

# Import objects to test
from flight_mech.wing import Wing, compute_chord_min_and_max_for_trapezoidal_wing
from flight_mech.airfoil import Airfoil

# Import test tools
from tests._common import check_value, output_folder

#############
# Constants #
#############


#########
# Tests #
#########

def test_create_3D_animation():
    wing = Wing()
    wing.y_array = np.linspace(0, 10, 100)
    wing.chord_length_array = np.linspace(3, 1, 100)
    wing.twisting_angle_array = np.linspace(0, 0.1, 100)
    wing.x_center_offset_array = np.zeros(100)
    airfoil = Airfoil("naca4412")
    wing.base_airfoil = airfoil
    wing.create_3D_animation(os.path.join(output_folder, "wing.gif"))
    wing.save_3D_shape(os.path.join(output_folder, "wing.stl"))


def test_plot_2D():
    wing = Wing()
    wing.y_array = np.linspace(0, 10, 100)
    wing.chord_length_array = np.linspace(3, 1, 100)
    wing.initialize()
    wing.plot_2D(hold_plot=True, save_path=os.path.join(
        output_folder, "wing.png"), clear_before_plot=True)

def test_surface():
    wing = Wing()
    wing.y_array = np.linspace(0, 10, 100)
    wing.chord_length_array = np.linspace(3, 1, 100)
    true_surface = 10 * (3 + 1) / 2
    assert wing.single_side_surface == true_surface

def test_compute_lift_and_induced_drag_coefficients_1():
    wing = Wing()
    nb_points_on_wing = 100
    airfoil = Airfoil("naca65210")
    true_aspect_ratio = 9
    chord_max = 72.6e-2  # m
    chord_min = 29e-2  # m
    true_wing_span = true_aspect_ratio * (chord_max + chord_min) / 2
    wing.y_array = np.linspace(0, true_wing_span / 2, nb_points_on_wing)
    wing.chord_length_array = np.linspace(
        chord_max, chord_min, nb_points_on_wing)
    wing.base_airfoil = airfoil
    wing.initialize()

    assert wing.aspect_ratio == true_aspect_ratio
    check_value(0.4, wing.taper_ratio)
    check_value(-1.2 * np.pi / 180, airfoil.compute_alpha_zero_lift())

    CL, CD = wing.compute_lift_and_induced_drag_coefficients(
        4 * np.pi / 180, 4)
    check_value(0.465, CL)

def test_compute_lift_and_induced_drag_coefficients_2():
    wing = Wing()
    nb_points_on_wing = 100
    airfoil = Airfoil("naca2412")
    true_aspect_ratio = 7.52
    true_surface = 16.3  # m2
    true_taper_ratio = 0.69

    chord_min, chord_max = compute_chord_min_and_max_for_trapezoidal_wing(
        true_surface, true_aspect_ratio, true_taper_ratio)
    wing_span = np.sqrt(true_aspect_ratio * true_surface)
    wing.y_array = np.linspace(0, wing_span / 2, nb_points_on_wing)
    wing.chord_length_array = np.linspace(
        chord_max, chord_min, nb_points_on_wing)
    wing.twisting_angle_array = np.linspace(
        1.5 * np.pi / 180, -1.5 * np.pi / 180, nb_points_on_wing)
    wing.base_airfoil = airfoil
    wing.initialize()

    check_value(true_aspect_ratio, wing.aspect_ratio)
    check_value(true_taper_ratio, wing.taper_ratio)
    check_value(true_surface, wing.reference_surface)
    check_value(-2.095 * np.pi / 180, airfoil.compute_alpha_zero_lift())

    CL, CD = wing.compute_lift_and_induced_drag_coefficients(
        0, 4)
    check_value(0.1593, CL)
    check_value(0.00166, CD)
