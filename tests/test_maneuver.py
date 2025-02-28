"""
Test file for the maneuver module.
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
import matplotlib.pyplot as plt

# Local imports #

# Import objects to test
from flight_mech.plane import Plane
from flight_mech.maneuver import TakeOffManeuver

# Import test tools
from tests._common import output_folder


#########
# Tests #
#########

def test_take_off():
    plane = Plane("su_27")

    # Adjust coefficient for take off configuration
    plane.C_D_0 = 0.023
    plane.C_L_alpha = 3.718
    plane.alpha_0 = -0.21 / plane.C_L_alpha
    plane.k = 0.1173
    plane.m_payload = 400 + 95 * 2 + 70  # kg

    # Check that the plane configuration is correct
    assert plane.m == 25556

    # Create the maneuver
    take_off_maneuver = TakeOffManeuver(
        plane_model=plane,
        rotation_speed=3 * np.pi / 180,
        rotation_sequence_trigger_speed=55.83,
        ground_friction_coefficient=0.03
    )

    assert take_off_maneuver.end_take_off_altitude == 30

    # Compute the evolution
    take_off_maneuver._nb_max_iterations = 100
    take_off_maneuver.dt = 0.25
    take_off_maneuver.compute_evolution()

    # Save the altitude graph
    take_off_maneuver.plot_graph("altitude", hold_plot=True, save_path=os.path.join(
        output_folder, "take_off_altitude.png"))

    assert np.isclose(take_off_maneuver.rotation_start_time, 7.25, rtol=0.1)
    assert np.isclose(take_off_maneuver.flight_start_time, 10.25, rtol=0.1)
    assert np.isclose(take_off_maneuver.x_coord_list[-1], 700, rtol=0.1)
