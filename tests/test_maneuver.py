"""
Test file for the maneuver module.
"""

###########
# Imports #
###########


# Python imports #

import os
import sys
sys.path.append(".")

# Dependencies #

import numpy as np
import matplotlib.pyplot as plt

# Local imports #

# Import objects to test
from flight_mech.maneuver import TakeOffManeuver

# Import test tools
from tests._common import output_folder


#########
# Tests #
#########
