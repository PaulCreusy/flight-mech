"""
Module to define common functions for the tests.
"""

###########
# Imports #
###########

import os

#############
# Constants #
#############

tolerance = 0.1
output_folder = os.path.join(os.path.dirname(__file__), "output")
data_folder = os.path.join(os.path.dirname(__file__), "data")

#############
# Functions #
#############

def check_value(true_value, test_value, tolerance=tolerance):
    assert abs(true_value - test_value) / true_value < tolerance
