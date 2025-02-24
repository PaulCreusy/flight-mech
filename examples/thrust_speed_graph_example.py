"""
Example file for the thrust-speed graph.
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

# Plot the thrust speed graph at several altitudes
plane.plot_gliding_TV_graph([0, 3000, 6000, 9000, 13000])
