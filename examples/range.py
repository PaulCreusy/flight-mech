"""
Example file for the range computations.
"""

##########
# Import #
##########

import sys
sys.path.append(".")
from flight_mech import Plane

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
