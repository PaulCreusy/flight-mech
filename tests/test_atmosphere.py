"""
Test file for the atmosphere module.
"""

###########
# Imports #
###########

# Python imports #

import sys
sys.path.append(".")

# Dependencies #

import numpy as np
import matplotlib.pyplot as plt

# Local imports #

# Import objects to test
from flight_mech.atmosphere import LinearAtmosphere, StandardAtmosphere

# Import test tools
from tests._common import check_value

#########
# Tests #
#########

def test_compute_density_from_altitude_linear():
    rho = LinearAtmosphere.compute_density_from_altitude(3000)
    check_value(0.9093, rho)

def test_compute_altitude_from_sigma_linear():
    altitude = LinearAtmosphere.compute_altitude_from_sigma(0.7891)
    check_value(2400, altitude)

def test_compute_sigma_from_altitude_linear():
    sigma = LinearAtmosphere.compute_sigma_from_altitude(4.2e3)
    check_value(0.6547, sigma)

def test_compute_density_from_altitude_standard():
    rho = StandardAtmosphere.compute_density_from_altitude(3000)
    check_value(0.9093, rho)

def test_compute_altitude_from_sigma_standard():
    altitude = StandardAtmosphere.compute_altitude_from_sigma(0.7891)
    check_value(2400, altitude)

def test_compute_sigma_from_altitude_standard():
    sigma = StandardAtmosphere.compute_sigma_from_altitude(4.2e3)
    check_value(0.6547, sigma)

def test_compute_temperature_from_altitude_standard():
    temperature = StandardAtmosphere.compute_temperature_from_altitude(9144)
    check_value(273 - 44.2, temperature)
