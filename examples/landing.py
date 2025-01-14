"""
Example file for the landing computations.
"""

##########
# Import #
##########

import sys
sys.path.append(".")
from src.flight_mech import Plane

###########
# Process #
###########

# Load the plane
plane = Plane("cessna_citation_III", "./plane_database")

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
