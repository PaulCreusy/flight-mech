"""
Example file for the range computations.
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

# Define the altitude
z = 6700  # m

# Print the results
print("max range at 6700m [m]", plane.compute_max_range_at_fixed_altitude(z))
speed = plane.compute_reference_speed(z)
print("v [m.s-1]", speed)
print("max range at fixed speed [m]",
      plane.compute_range_at_fixed_speed(speed, f=plane.f_max))
