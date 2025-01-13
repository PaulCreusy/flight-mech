"""
Example file for the take-off computations.
"""

##########
# Import #
##########

import sys
sys.path.append(".")
from src.plane import Plane

###########
# Process #
###########

# Load the plane
plane = Plane("cessna_citation_III", "./plane_database")

# Define the friction coefficient of the runway
mu = 0.02

# Define the take off altitude
z = 0

# Print the results
print("ground effect coef", plane.ground_effect_coefficient)
print("take off speed [m.s-1]", plane.compute_take_off_speed(z))
print("take off distance [m]",
      plane.compute_take_off_distance_with_friction(z, mu))
