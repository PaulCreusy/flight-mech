"""
Test file for the aerodynamics module.
"""

###########
# Imports #
###########

# Python imports #

import sys
sys.path.append(".")

# Dependencies #

import numpy as np

# Local imports #

# Import objects to test
from flight_mech.aerodynamics import (
    compute_blasius_rectangle_drag,
    compute_polhausen_linear_drag,
    compute_blasius_linear_drag
)


#########
# Tests #
#########

def test_compute_rectangle_drag_blasius():
    rectangle_drag = compute_blasius_rectangle_drag(
        width=3, length=1, rho=1.23, nu=1.46e-5, velocity=2)
    assert np.isclose(rectangle_drag, 2.65e-2, rtol=0.001).all()

def test_compute_polhausen_linear_drag():
    nb_points = 3000
    x_array = np.linspace(0, 1, nb_points)
    velocity_array = 2 * np.ones(nb_points)
    rho = 1.23  # kg.m-3
    nu = 1.46e-5  # m2.s-1
    blasius_linear_drag = compute_blasius_linear_drag(
        x_array, velocity_array, rho, nu)
    polhausen_linear_drag = compute_polhausen_linear_drag(
        x_array, velocity_array, rho, nu)
    assert np.isclose((0.343 / 0.332) * blasius_linear_drag,
                      polhausen_linear_drag, rtol=0.1).all()
