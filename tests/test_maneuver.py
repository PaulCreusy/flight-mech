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
from flight_mech.maneuver import TakeOffManeuver, LandingManeuver

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
    plane.k = 0.11732
    plane.m_payload = 400 + 95 * 2 + 70  # kg
    plane.update_P(force=True)

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
    take_off_maneuver.nb_max_iterations = 100
    take_off_maneuver.dt = 0.25
    take_off_maneuver.compute_evolution()

    # Save the altitude graph
    take_off_maneuver.plot_graph("altitude", hold_plot=True, save_path=os.path.join(
        output_folder, "take_off_altitude.png"))

    assert np.isclose(take_off_maneuver.rotation_start_time, 7.25, rtol=0.1)
    assert np.isclose(take_off_maneuver.flight_start_time, 10.25, rtol=0.1)
    assert np.isclose(
        take_off_maneuver.ground_distance_list[-1], 700, rtol=0.1)

def test_landing():
    plane = Plane("su_27")

    # Adjust coefficient for take off configuration
    plane.C_D_0 = 0.0394
    plane.C_L_alpha = 3.718
    plane.alpha_0 = -0.428 / plane.C_L_alpha
    plane.k = 0.11732
    plane.m_payload = 400 + 95 * 2 + 70  # kg
    plane.m_fuel = plane.m_fuel * 0.1
    plane.update_P(force=True)

    # Check that the plane configuration is correct
    assert np.isclose(plane.m, 18014.8, rtol=0.001)

    # Create the maneuver
    landing_maneuver = LandingManeuver(
        plane_model=plane,
        rotation_speed=4 * np.pi / 180,
        initial_velocity=62.34,
        parachute_drag_coefficient=0.6,
        parachute_reference_surface=15,
        ground_friction_coefficient=0.03,
        braking_friction_coefficient=0.5
    )

    # Compute the evolution
    landing_maneuver.nb_max_iterations = 800
    landing_maneuver.dt = 0.1
    landing_maneuver.compute_evolution()

    landing_maneuver.plot_graph("velocity", use_grid=True, hold_plot=True, save_path=os.path.join(
        output_folder, "landing_velocity.png"))

    assert np.isclose(landing_maneuver.braking_start_time, 3, rtol=0.05)
    assert np.isclose(landing_maneuver.time_list[-1], 12.9, rtol=0.1)
    assert np.isclose(
        landing_maneuver.ground_distance_list[-1], 434, rtol=0.1)
