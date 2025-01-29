"""
Example file for the take-off computations.
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

# Define the friction coefficient of the runway
mu = 0.02

# Define the take off altitude
z = 0

# Print the results
print("ground effect coef", plane.ground_effect_coefficient)
print("take off speed [m.s-1]", plane.compute_take_off_speed(z))
print("take off distance [m]",
      plane.compute_take_off_distance_with_friction(z, mu))
