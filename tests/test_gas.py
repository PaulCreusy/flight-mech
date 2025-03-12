"""
Test file for the gas module.
"""

###########
# Imports #
###########

# Python imports #

import sys
sys.path.append(".")

# Dependencies #

import numpy as np
import matplotlib.pyplot as plt

# Local imports #

# Import objects to test
from flight_mech.gas import MonoatomicPerfectGas, DiatomicPerfectGas, GasMixture


#########
# Tests #
#########

def test_monoatomic_perfect_gas():
    helium = MonoatomicPerfectGas("He")
    helium.pressure = 1e5  # Pa
    helium.temperature = 273  # K

    assert np.isclose(helium.density, 0.1786, rtol=0.05)
    assert np.isclose(helium.sound_velocity, 972, rtol=0.05)
    assert np.isclose(helium.dynamic_viscosity, 1.87e-5, rtol=0.05)


def test_diatomic_perfect_gas():
    hydrogen = DiatomicPerfectGas("H2", pressure=1e5, temperature=293)

    assert np.isclose(hydrogen.density, 0.08988, rtol=0.1)
    assert np.isclose(hydrogen.Cp, 14266, rtol=0.05)

    hydrogen.T = 273 + 27
    assert np.isclose(hydrogen.sound_velocity, 1310, rtol=0.05)

def test_gas_mixture():
    oxygen = DiatomicPerfectGas("O2")
    nitrogen = DiatomicPerfectGas("N2")
    gas_model_dict = {
        "O2": oxygen,
        "N2": nitrogen
    }
    gas_molar_fraction_dict = {
        "O2": 0.21,
        "N2": 0.79
    }

    air = GasMixture(gas_model_dict, gas_molar_fraction_dict)
    air.temperature = 283
    air.P = 1e5

    assert np.isclose(air.density, 1.225, rtol=0.01)
    assert np.isclose(air.sound_velocity, 330, rtol=0.05)
    assert np.isclose(air.Cp, 1005, rtol=0.01)
    assert np.isclose(air.dynamic_viscosity, 17.26e-6, rtol=0.05)

    air.temperature = 750
    assert np.isclose(air.density, 0.471, rtol=0.05)
