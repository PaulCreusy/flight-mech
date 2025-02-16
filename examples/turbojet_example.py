"""
Example file for the turbojet class.
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

from flight_mech.turbine import TurbojetSingleBody

#################
# Configuration #
#################

# Define the flight characteristics
cruise_altitude = 11800  # m
cruise_mach = 0.78
max_mach = 0.82
max_altitude = 12500  # m


# Define the turbojet
turbojet = TurbojetSingleBody()

# Set the design inputs
turbojet.compressor_efficiency = 0.86
turbojet.turbine_efficiency = 0.9
turbojet.OPR_design = 10
turbojet.T4_max = 1700  # K

# Set the operating conditions
turbojet.M0 = 0
turbojet.ambient_pressure = 101325  # Pa
turbojet.ambient_temperature = 285  # K

##########
# Design #
##########

# Tune A4* to obtain 7500N of thrust
turbojet.tune_A4_star_for_desired_thrust(7500)
print("A4* at 7500N [m2]:", turbojet.A4_star)

#############
# Operation #
#############

# Switch to operation mode
turbojet.mode = "operation"

# Set the number of points for the plots
nb_points = 100

# T4 influence #

# Compute thrust evolution with T4 instruction
T4_array = np.linspace(0, turbojet.T4_max, nb_points)
thrust_array = np.zeros(nb_points)
for i, T4 in enumerate(T4_array):
    turbojet.T4_instruction = T4
    turbojet.tune_current_OPR()
    thrust_array[i] = turbojet.thrust

# Plot thrust evolution with velocity
plt.plot(T4_array, thrust_array)
plt.xlabel("T4 [K]")
plt.ylabel("Thrust [N]")
plt.title("Thrust evolution with T4")
plt.show()

# Cruise #

# Switch to cruise conditions
# T4 not realistic but to see what happens at full power
turbojet.T4_instruction = turbojet.T4_max
turbojet.altitude = cruise_altitude
turbojet.M0 = cruise_mach

# Compute thrust in cruise
turbojet.tune_current_OPR()
print("cruise thrust [N]:", turbojet.thrust)

# Velocity influence #

# Compute thrust evolution with velocity
mach_array = np.linspace(0, max_mach, nb_points)
thrust_array = np.zeros(nb_points)
for i, mach in enumerate(mach_array):
    turbojet.M0 = mach
    turbojet.tune_current_OPR()
    thrust_array[i] = turbojet.thrust

# Plot thrust evolution with velocity
plt.plot(mach_array, thrust_array)
plt.xlabel("M0")
plt.ylabel("Thrust [N]")
plt.title("Thrust evolution with Mach")
plt.show()

# Altitude influence #

# Compute thrust evolution with altitude
turbojet.M0 = cruise_mach
altitude_array = np.linspace(0, max_altitude, nb_points)
thrust_array = np.zeros(nb_points)
for i, altitude in enumerate(altitude_array):
    turbojet.altitude = altitude
    turbojet.tune_current_OPR()
    thrust_array[i] = turbojet.thrust

# Plot thrust evolution with altitude
plt.plot(altitude_array, thrust_array)
plt.xlabel("Altitude [m]")
plt.ylabel("Thrust [N]")
plt.title("Thrust evolution with altitude")
plt.show()
