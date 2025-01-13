"""
Example file for the thrust-speed graph.
"""

##########
# Import #
##########

import sys
import numpy as np
sys.path.append(".")
from src.plane import Plane

###########
# Process #
###########

# Load the plane
plane = Plane("cessna_citation_III", "./plane_database")

# Set the missing lift coefficient
plane.a = 2 * np.pi

# Plot the thrust speed graph at several altitudes
plane.plot_gliding_TV_graph([0, 3000, 6000, 9000, 13000])
