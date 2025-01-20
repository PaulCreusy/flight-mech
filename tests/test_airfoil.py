"""
Test file for the airfoil module.
"""

###########
# Imports #
###########

# Python imports #

import os
import sys
sys.path.append(".")

# Local imports #

# Import objects to test
from flight_mech.airfoil import Airfoil

#############
# Constants #
#############

airfoil_database_path = "./flight_mech/airfoil_database/"

#########
# Tests #
#########

def test_load_selig_file():
    airfoil = Airfoil()
    airfoil.load_selig_file(os.path.join(
        airfoil_database_path, "FX 62-K-153.txt"), 1)


test_load_selig_file()
