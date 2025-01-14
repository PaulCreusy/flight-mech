"""
Example file for the polar graph.
"""

##########
# Import #
##########

import sys
sys.path.append(".")
import numpy as np
from src.flight_mech import Plane

###########
# Process #
###########

# Load the plane
plane = Plane("cessna_citation_III", "./plane_database")

# Set the missing lift coefficient
plane.a = 2 * np.pi

# Plot the polar graph
plane.plot_polar_graph()
