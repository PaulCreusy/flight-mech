"""
Module to define fuel characteristics.
"""

###########
# Imports #
###########

# Python imports #

from enum import Enum


###########
# Classes #
###########

class FuelType:
    density: float  # kg.m-3
    lower_heating_value: float  # J

    def __init__(self, density: float, lower_heating_value: float):
        self.density = density
        self.lower_heating_value = lower_heating_value

class FuelTypes(FuelType, Enum):
    METHANE = 0.670, 50.1e6
    DIESEL = 0.91e3, 43e6
    H2 = 0.088, 120.1e6
    GASOLINE = 0.75e3, 41e6
    KEROSENE = 0.8e3, 44e6
    JET_A = 807.5, 42.8e6
    JP_4 = 776.5, 42.8e6
    JP_5 = 816.5, 42.8e6
