"""
Test file for the airfoil module.
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

# Local imports #

# Import objects to test
from flight_mech.airfoil import Airfoil

# Import test tools
from tests._common import check_value

#############
# Constants #
#############

airfoil_database_path = "./flight_mech/airfoil_database/"

#########
# Tests #
#########

def test_load_selig_file():
    airfoil = Airfoil()
    airfoil.load_selig_file(os.path.join(
        airfoil_database_path, "fx62k153.txt"), 1)

def test_import_from_airfoiltools():
    airfoil = Airfoil()
    airfoil.import_from_airfoiltools(
        "NACA", max_thickness=16, min_thickness=15, maximise_glide_ratio_at_reynolds="50k")

def test_re_interpolate():
    airfoil = Airfoil()
    airfoil.load_database_airfoil("fx62k153")
    airfoil.re_interpolate(np.linspace(0, 1, 1000))

def test_max_thickness():
    airfoil = Airfoil()
    airfoil.load_database_airfoil("naca4412")
    check_value(0.12, airfoil.max_thickness)
    check_value(0.3, airfoil.max_thickness_location)

def test_chord_length():
    airfoil = Airfoil()
    airfoil.load_database_airfoil("fx62k153")
    thickness_1 = airfoil.max_thickness
    airfoil.chord_length = 2
    assert np.max(airfoil.x_array) == 2
    assert airfoil.max_thickness == thickness_1

def test_max_camber():
    airfoil = Airfoil("fx62k153")
    check_value(0.041, airfoil.max_camber)
    check_value(0.629, airfoil.max_camber_location)
