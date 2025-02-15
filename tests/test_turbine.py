"""
Test file for the turbine module.
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
import pandas as pd
import matplotlib.pyplot as plt

# Local imports #

# Import objects to test
from flight_mech.turbine import TurbojetSingleBody, VARIABLE_TO_CODE

# Import test tools
from tests._common import check_value, output_folder


#########
# Tests #
#########

def test_compute_thrust():
    turbojet = TurbojetSingleBody()
    turbojet.A4_star = 5e-2
    assert np.isclose(turbojet.thrust, 46117.878152741425, rtol=0.001)

def test_tune_A4_star_for_desired_thrust():
    turbojet = TurbojetSingleBody()
    turbojet.tune_A4_star_for_desired_thrust(7500)
    assert np.isclose(turbojet.A4_star, 0.008131337010233023, rtol=0.001)

def test_plot_graph():
    turbojet = TurbojetSingleBody()
    for variable in VARIABLE_TO_CODE:
        turbojet.plot_graph(variable, hold_plot=True, clear_before_plot=True, save_path=os.path.join(
            output_folder, f"{variable}_graph.png"))
