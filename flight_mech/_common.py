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
# Constants #
#############

LAZY_MEMORY_LIMIT = 10

#############
# Functions #
#############

def plot_graph(
        x_array: np.ndarray,
        y_array: np.ndarray,
        data_label: str | None = None,
        title: str | None = None,
        use_grid: bool = False,
        use_legend: bool = False,
        save_path: str | None = None,
        hold_plot: bool = False,
        clear_before_plot: bool = False,
        axis_type: str | None = None,
        x_label: str | None = None,
        y_label: str | None = None):
    """
    Plot a graph with the selected arguments.

    Parameters
    ----------
    x_array : np.ndarray
        X data.
    y_array : np.ndarray
        Y data.
    data_label : str | None, optional
        Data label, by default None
    title : str | None, optional
        Title, by default None
    use_grid : bool, optional
        Indicate wether to use grid, by default False
    save_path : str | None, optional
        Path to save the figure, by default None
    hold_plot : bool, optional
        Indicate wether the graph is kept or plot, by default False
    clear_before_plot : bool, optional
        Indicate wether the graph is cleared before plotting, by default False
    axis_type : str | None, optional
        Type of axis, by default None
    x_label : str | None, optional
        X data label, by default None
    y_label : str | None, optional
        Y data label, by default None
    """

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

    # Enable legend if needed
    if use_legend:
        plt.legend()

    # Save figure if needed
    if save_path is not None:
        plt.savefig(save_path)

    # Show it if needed
    if not hold_plot:
        plt.show()

def lazy(func):
    memory = {}

    def inner(arg):
        current_hash = hash(arg)
        if current_hash in memory:
            res = memory[current_hash]
            return res
        if len(memory) > LAZY_MEMORY_LIMIT:
            memory.clear()
        res = func(arg)
        memory[current_hash] = res
        return res

    return inner

def hash_dict(dict_to_hash: dict):
    res = 0
    for key in dict_to_hash:
        res += (hash(key) * hash(dict_to_hash[key])) % int(1e12)
    return res
