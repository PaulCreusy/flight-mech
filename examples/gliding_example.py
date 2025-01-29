"""
Example file for the gliding computations.
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
z = 300

# Print the results
print("R max [m]", plane.compute_max_gliding_range(z))
print("gamma min [°]", plane.compute_min_descent_gliding_slope() * 180 / np.pi)
