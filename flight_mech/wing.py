"""
Module to analyse wing aerodynamic properties.
"""

###########
# Imports #
###########

# Python imports #

import os
from typing import Literal

# Dependencies #

import numpy as np
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt
from copy import deepcopy
try:
    import pyvista as pv
except:
    raise Warning(
        "Pyvista is not detected, please install it before using the 3D visualization functions.")

# Local imports #

from flight_mech._common import plot_graph
from flight_mech.airfoil import Airfoil

#############
# Constants #
#############

# Define a default plane database location
default_wing_database = os.path.join(
    os.path.dirname(__file__), "wing_database")

#############
# Functions #
#############

def check_pyvista_import():
    try:
        pv.__version__
    except:
        raise ImportError(
            "You need to install pyvista before using 3D visualization functions. You can do it with: 'pip install pyvista'")

###########
# Classes #
###########

class Wing:
    name: str | None = None
    _y_array: np.ndarray | None = None
    chord_length_array: np.ndarray | None = None
    twisting_angle_array: np.ndarray | None = None
    x_center_offset_array: np.ndarray | None = None
    base_airfoil: Airfoil | None = None

    @property
    def y_array(self) -> np.ndarray:
        """
        Array of y coordinates used for the wing definition.
        """
        return self._y_array

    @y_array.setter
    def y_array(self, value: np.ndarray):
        if self.twisting_angle_array is not None and self.chord_length_array is not None and self.x_center_offset_array is not None:
            if self._y_array.size == self.chord_length_array.size:
                self.re_interpolate(value, update_y_array=False)
            self._y_array = value
        else:
            self._y_array = value
        self._chord_length = np.max(value)

    @property
    def single_side_surface(self):
        """
        Surface of a single wing.
        """
        surface = np.trapezoid(self.chord_length_array, self.y_array)
        return surface

    @property
    def reference_surface(self):
        """
        Reference surface of the wings of the plane. Equal to 2 times the surface of a single wing.
        """
        reference_surface = self.single_side_surface * 2
        return reference_surface

    @property
    def wing_span(self):
        """
        Wing span.
        """
        wing_span = np.max(self.chord_length_array)
        return wing_span

    @property
    def aspect_ratio(self):
        """
        Aspect ratio.
        """
        aspect_ratio = pow(self.wing_span, 2) / self.reference_surface
        return aspect_ratio

    def initialize(self):
        if self.y_array is None:
            raise ValueError(
                "Unable to initialize the wing, the y array is not defined.")
        if self.twisting_angle_array is None:
            self.twisting_angle_array = np.zeros(self.y_array.shape)
        if self.x_center_offset_array is None:
            self.x_center_offset_array = np.zeros(self.y_array.shape)
        if self.base_airfoil is None:
            self.base_airfoil = Airfoil("naca4412")

    def check_initialization(self):
        if self.y_array is None:
            raise ValueError(
                "The y array is not defined, unable to proceed. Please define it first.")
        if self.chord_length_array is None:
            raise ValueError(
                "The chord length array is not defined, unable to proceed. Please define it first.")
        if self.twisting_angle_array is None:
            raise ValueError(
                "The twisting angle array is not defined, unable to proceed. Please initialize the wing with wing.initialize() or define it first.")
        if self.x_center_offset_array is None:
            raise ValueError(
                "The x center offset array is not defined, unable to proceed. Please initialize the wing with wing.initialize() or define it first.")
        if self.base_airfoil is None:
            raise ValueError(
                "The base airfoil is not defined, unable to proceed. Please initialize the wing with wing.initialize() or define it first.")

    def re_interpolate(self, new_y_array: np.ndarray, update_y_array: bool = True):
        pass

    def plot_2D(self, save_path: str | None = None, hold_plot: bool = False, clear_before_plot: bool = False):
        """
        Plot the shape of the wing in 2D.

        Parameters
        ----------
        save_path : str | None, optional
            Path to save the figure, by default None
        hold_plot : bool, optional
            Indicate wether to keep or plot the figure, by default False
        clear_before_plot : bool, optional
            Indicate wether to clear or not the plot before, by default False
        """

        # Allocate arrays to define the leading edge and trailing edge of the wing
        leading_edge_x_array = self.x_center_offset_array + (
            self.chord_length_array[0] - self.chord_length_array) / 2
        trailing_edge_x_array = self.chord_length_array + self.x_center_offset_array + (
            self.chord_length_array[0] - self.chord_length_array) / 2

        # Group to create a contour
        x_contour_array = np.concatenate(
            (leading_edge_x_array, trailing_edge_x_array[::-1], [leading_edge_x_array[0]]), axis=0)
        y_contour_array = np.concatenate(
            (self.y_array, self.y_array[::-1], [self.y_array[0]]), axis=0)

        # Plot the graph
        plot_graph(
            y_contour_array,
            x_contour_array,
            data_label=self.name,
            title="Wing visualization",
            use_grid=True,
            save_path=save_path,
            hold_plot=hold_plot,
            clear_before_plot=clear_before_plot,
            x_label="y",
            y_label="x",
            axis_type="equal"
        )

    def create_3D_animation(self, output_path: str, nb_points_airfoil=50, nb_frames=60, time_step=0.05):
        """
        Create a rotating 3D animation of the wing.

        Parameters
        ----------
        output_path : str
            Path of the output gif.
        nb_points_airfoil : int, optional
            Number of points to use for the airfoil, by default 50
        nb_frames : int, optional
            Number of frames in the animation, by default 60
        time_step : float, optional
            Time step between each frame, by default 0.05
        """

        # Check if pyvista is imported
        check_pyvista_import()

        # Create the wing surface
        wing_surface = self.create_wing_3D_surface(nb_points_airfoil)

        pl = pv.Plotter(off_screen=True)
        pl.add_mesh(wing_surface)
        pl.show(auto_close=False)
        path = pl.generate_orbital_path(
            n_points=nb_frames, shift=wing_surface.length, factor=3.0)
        if not output_path.endswith(".gif"):
            output_path = output_path + ".gif"
        pl.open_gif(output_path)
        pl.orbit_on_path(path, write_frames=True, step=time_step)
        pl.close()

    def create_wing_3D_surface(self, nb_points_airfoil=50):
        """
        Create a wing surface PyVista object for 3D plotting.

        Parameters
        ----------
        nb_points_airfoil : int, optional
            Number of points to use for the airfoil, by default 50

        Returns
        -------
        pv.PointSet
            Pyvista 3D surface.
        """

        # Check if pyvista is imported
        check_pyvista_import()

        # Create a display airfoil with normalized chord and less points
        display_airfoil = deepcopy(self.base_airfoil)
        display_airfoil.chord_length = 1
        display_airfoil.re_interpolate_with_cosine_distribution(
            nb_points_airfoil // 2)

        # Create an array containing all the points
        points_array = np.zeros((nb_points_airfoil * self.y_array.size * 2, 3))
        for i in range(self.y_array.size):
            ratio = self.chord_length_array[i] / display_airfoil.chord_length
            airfoil_x_array, airfoil_z_array = display_airfoil.get_rotated_selig_arrays(
                self.twisting_angle_array[i])
            points_index = i * nb_points_airfoil
            points_array[points_index:points_index + nb_points_airfoil, 0] =\
                airfoil_x_array[:-1] * ratio + \
                self.x_center_offset_array[i] + (
                    self.chord_length_array[0] - self.chord_length_array[i]) / 2
            points_array[points_index:points_index +
                         nb_points_airfoil, 1] = self.y_array[i]
            points_array[points_index:points_index + nb_points_airfoil,
                         2] = airfoil_z_array[:-1] * ratio

        # Create a mesh from the points
        point_cloud = pv.PolyData(points_array)
        wing_surface = point_cloud.delaunay_3d()

        return wing_surface

    def plot_3D(self, nb_points_airfoil=50, show_symmetric_wing: bool = False):
        """
        Plot the shape of the wing in 3D.
        """

        # Create the 3D surface
        wing_surface = self.create_wing_3D_surface(nb_points_airfoil)

        # Create the plotter object
        p = pv.Plotter()

        # Create a symmetry if needed
        if show_symmetric_wing:
            wing_reflected = wing_surface.reflect((0, 1, 0))
            p.add_mesh(wing_reflected)

        # Plot
        p.add_mesh(wing_surface)
        p.show()

    def compute_lift_and_drag(alpha: float, nb_points_fourrier: int):

        # Create theta positions to interpolate
        theta_interpolation_pos = np.linspace(
            np.pi / 2 / nb_points_fourrier, np.pi / 2, nb_points_fourrier)

        # Interpolate the wing values on the theta positions
        y_on_theta_array = ...
