"""
Test file for the atmosphere module.
"""

###########
# Imports #
###########

# Python imports #

import sys
sys.path.append(".")

# Local imports #

# Import objects to test
from flight_mech.atmosphere import compute_air_density_from_altitude, compute_altitude_from_sigma, compute_sigma_from_altitude

# Import test tools
from tests._common import check_value

#########
# Tests #
#########

def test_compute_air_density_from_altitude():
    rho = compute_air_density_from_altitude(3000)
    check_value(0.9093, rho)

def test_compute_altitude_from_sigma():
    altitude = compute_altitude_from_sigma(0.7891)
    check_value(2400, altitude)

def test_compute_sigma_from_altitude():
    sigma = compute_sigma_from_altitude(4.2e3)
    check_value(0.6547, sigma)
