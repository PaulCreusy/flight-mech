"""
Example file for the ascension computations.
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
plane = Plane("cessna_172")

# Compute the fmax and CL at fmax
print("C_L_f_max", plane.C_L_f_max)
print("fmax", plane.f_max)

# Compute the speed interval at 8000 meters
plane.m_fuel = 136.26  # kg
plane.update_variables(True)
print("reference speed at 8000m [m.s-1]", plane.compute_reference_speed(8000))
print("speed interval at 8000m [m.s-1]",
      plane.compute_velocity_interval_for_fixed_thrust(8000))
print("stall speed at 8000m [m.s-1]",
      plane.compute_stall_speed(8000, C_L_max=1.5))

# Compute the ascension speed and slope at sea level
plane.m_fuel = 0  # kg
plane.update_variables(True)
print("max ascension speed [m.s-1]", plane.compute_max_ascension_speed(z=0))
print("reference speed at 0m [m.s-1]", plane.compute_reference_speed(z=0))
print("max slope at 0m [%]", plane.compute_max_ascension_slope(z=0))
