"""
Module to define common functions for all flight-mech modules.
"""

###########
# Imports #
###########

# Dependencies #

import numpy as np
import matplotlib.pyplot as plt

#############
# Functions #
#############

def plot_graph(
        x_array: np.ndarray,
        y_array: np.ndarray,
        data_label: str | None = None,
        title: str | None = None,
        use_grid: bool = False,
        save_path: str | None = None,
        hold_plot: bool = False,
        clear_before_plot: bool = False,
        axis_type: str | None = None,
        x_label: str | None = None,
        y_label: str | None = None):

    # Clear plot if needed
    if clear_before_plot:
        plt.cla()
        plt.clf()

    # Add data
    plt.plot(x_array, y_array, label=data_label)

    # Add labels
    if x_label is not None:
        plt.xlabel(x_label)
    if y_label is not None:
        plt.ylabel(y_label)

    # Add title
    if title is not None:
        plt.title(title)

    # Set axis type
    if axis_type is not None:
        plt.axis(axis_type)

    # Enable grid if needed
    if use_grid:
        plt.grid()

    # Save figure if needed
    if save_path is not None:
        plt.savefig(save_path)

    # Show it if needed
    if not hold_plot:
        plt.show()
