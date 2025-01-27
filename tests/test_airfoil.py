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
import pandas as pd
import matplotlib.pyplot as plt

# Local imports #

# Import objects to test
from flight_mech.airfoil import Airfoil

# Import test tools
from tests._common import check_value

#############
# Constants #
#############

airfoil_database_path = "./flight_mech/airfoil_database/"
airfoil_tools_predictions_naca4412 = pd.read_csv(
    "tests/data/xf-naca4412-il-500000.csv", skiprows=10)
coefficient_error_threshold = 0.05

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
    check_value(thickness_1, airfoil.max_thickness)

def test_max_camber():
    airfoil = Airfoil("fx62k153")
    check_value(0.041, airfoil.max_camber)
    check_value(0.629, airfoil.max_camber_location)

def test_compute_airfoil_fourrier_coefficients():
    airfoil = Airfoil("n11h9")
    airfoil.compute_airfoil_fourrier_coefficients()
    check_value(0.1049, airfoil._a0, tolerance=0.2)
    check_value(0.2501, airfoil._a1, tolerance=0.2)
    check_value(0.2388, airfoil._a2, tolerance=0.2)

def test_compute_lift_coefficient():
    airfoil = Airfoil("naca4412")
    alpha_deg = airfoil_tools_predictions_naca4412["Alpha"]
    mask = np.abs(alpha_deg) <= 7.5
    alpha = alpha_deg[mask] * np.pi / 180
    CL_airfoil_tools = airfoil_tools_predictions_naca4412["Cl"][mask]
    CL = airfoil.compute_lift_coefficient(alpha)
    max_diff = np.max(np.abs(CL - CL_airfoil_tools))
    assert max_diff < coefficient_error_threshold

def test_compute_momentum_coefficient():
    airfoil = Airfoil("naca4412")
    alpha_deg = airfoil_tools_predictions_naca4412["Alpha"]
    mask = np.abs(alpha_deg) <= 7.5
    alpha = alpha_deg[mask] * np.pi / 180
    Cm_airfoil_tools = airfoil_tools_predictions_naca4412["Cm"][mask]
    Cm = airfoil.compute_momentum_coefficient_at_aero_center(alpha)
    max_diff = np.max(np.abs(Cm - Cm_airfoil_tools))
    assert max_diff < coefficient_error_threshold
