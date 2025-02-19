"""
Test file for the motor module.
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
from flight_mech.motor import ElectricalMotor

# Import test tools
from tests._common import output_folder


#########
# Tests #
#########


def test_compute_electromotive_force():
    motor = ElectricalMotor("electron_60")
    motor.U = 240  # V
    motor.I = 34  # A
    assert np.isclose(motor.electromotive_force, 229.6347, rtol=0.001)

def test_compute_rotation_speed():
    motor = ElectricalMotor("electron_60")
    motor.U = 240  # V
    motor.I = 34  # A
    assert np.isclose(motor.rotation_speed * 60 /
                      (2 * np.pi), 3214.8858, rtol=0.001)

def test_compute_efficiency():
    motor = ElectricalMotor("electron_60")
    motor.U = 240  # V
    motor.I = 34  # A
    assert np.isclose(motor.efficiency, 0.815883676, rtol=0.001)

def test_compute_power():
    motor = ElectricalMotor("electron_60")
    motor.U = 240  # V
    motor.I = 34  # A
    assert np.isclose(motor.power, 6657.6108, rtol=0.001)

def test_compute_torque():
    motor = ElectricalMotor("electron_60")
    motor.U = 240  # V
    motor.I = 34  # A
    assert np.isclose(motor.torque, 19.77535254, rtol=0.001)

def test_plot_graph():
    motor = ElectricalMotor("electron_60")
    motor.plot_graph("efficiency", nb_points=1000, hold_plot=True, save_path=os.path.join(
        output_folder, "motor_efficiency.png"))

def test_compute_I_from_rotation_speed():
    motor = ElectricalMotor("electron_60")
    motor.U = 240  # V
    motor.rotation_speed = 2726.9886 * (2 * np.pi) / 60
    assert np.isclose(motor.I, 148, rtol=0.001)

def test_compute_I_from_power():
    motor = ElectricalMotor("electron_60")
    motor.U = 240  # V
    motor.power = 18058.83
    assert np.isclose(motor.I, 90, rtol=0.001)

def test_compute_I_from_torque():
    motor = ElectricalMotor("electron_60")
    motor.U = 240  # V
    motor.torque = 90.69168952
    assert np.isclose(motor.I, 138, rtol=0.001)

def test_compute_I_at_max_efficiency():
    motor = ElectricalMotor("electron_60")
    motor.U = 240  # V
    I_at_max_efficiency = motor.compute_I_at_max_efficiency()
    assert np.isclose(I_at_max_efficiency, 62, rtol=0.03)
