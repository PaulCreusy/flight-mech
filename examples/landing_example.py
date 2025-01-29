"""
Example file for the landing computations.
"""

##########
# Import #
##########

import os
import sys
sys.path.append(".")
import numpy as np
import matplotlib.pyplot as plt

# If running inside pytest subprocess, mock plt.show
if os.getenv("PYTEST_RUNNING"):
    plt.show = lambda: None

from flight_mech.plane import Plane

###########
# Process #
###########

# Load the plane
plane = Plane("cessna_citation_III")

# Define the friction coefficient of the landing strip
mu = 0.1

# Define the landing altitude
z = 0

# Update some coefficients to set the landing configuration
plane.C_D_0 = plane.C_D_0 * 1.1
plane.m_fuel = 0
plane.C_L_max = 2.5
plane.update_variables(force=True)

# Print the results
print("weight [N]", plane.P)
print("landing speed [m.s-1]", plane.compute_landing_speed(z))
print("landing distance [m]", plane.compute_landing_distance(z, mu, C_L=0))
