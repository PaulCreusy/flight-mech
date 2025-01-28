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
from flight_mech.wing import Wing
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
