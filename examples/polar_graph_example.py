"""
Example file for the polar graph.
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

# Set the missing lift coefficient
plane.C_L_alpha = 2 * np.pi

# Plot the polar graph
plane.plot_polar_graph()
