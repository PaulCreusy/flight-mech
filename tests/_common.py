"""
Module to define common functions for the tests.
"""

#############
# Constants #
#############

tolerance = 0.1

#############
# Functions #
#############

def check_value(true_value, test_value, tolerance=tolerance):
    assert abs(true_value - test_value) / true_value < tolerance
