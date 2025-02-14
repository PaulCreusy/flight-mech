"""
Module to analyse wing aerodynamic properties.
"""

###########
# Imports #
###########

# Python imports #

import os
from typing import Literal
from copy import deepcopy

# Dependencies #

import numpy as np
from scipy.interpolate import make_interp_spline

try:
    import pyvista as pv
except Exception as exc:
    raise Warning(
        "Pyvista is not detected, please install it before using the 3D visualization functions.") from exc

# Local imports #

from flight_mech._common import plot_graph
from flight_mech.aerodynamics import compute_linear_drag
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
    """
    Verify that pyvista is successfully imported.

    Raises
    ------
    ImportError
        Raise error if pyvista cannot be accessed.
    """

    try:
        pv.__version__
    except Exception as exc:
        raise ImportError(
            "You need to install pyvista before using 3D visualization functions. You can do it with: 'pip install pyvista'") from exc

def convert_y_to_theta(y_array: np.ndarray, wing_span: float):
    """
    Convert a y array to theta.

    Parameters
    ----------
    y_array : np.ndarray
        Values of y coordinates.
    wing_span : float
        Wing span.

    Returns
    -------
    np.ndarray
        Theta array.
    """

    theta_array = np.acos(2 * y_array / wing_span)

    return theta_array

def convert_theta_to_y(theta_array: np.ndarray, wing_span: float):
    """
    Convert a theta array to y.

    Parameters
    ----------
    theta_array : np.ndarray
        Values of theta coordinates.
    wing_span : float
        Wing span.

    Returns
    -------
    np.ndarray
        Y coordinates array.
    """

    y_array = (wing_span / 2) * np.cos(theta_array)

    return y_array

def compute_chord_min_and_max_for_trapezoidal_wing(reference_surface: float, aspect_ratio: float, taper_ratio: float):
    """
    Compute the chord min and max values for a trapezoidal wing.

    Parameters
    ----------
    reference_surface : float
        Reference surface of the wing.
    aspect_ratio : float
        Aspect ratio of the wing.
    taper_ratio : float
        Taper ratio of the wing.

    Returns
    -------
    tuple[float,float]
        Tuple containing the chord min and max values.
    """

    wing_span = np.sqrt(aspect_ratio * reference_surface)
    max_chord = 2 * reference_surface / (wing_span * (1 + taper_ratio))
    min_chord = max_chord * taper_ratio

    return min_chord, max_chord

###########
# Classes #
###########

class Wing:
    """
    Class to define a wing and compute its characteristics.
    """

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

    @property
    def leading_edge_x_array(self) -> np.ndarray:
        """
        Array of x coordinates of the leading edge.
        """
        leading_edge_x_array = self.x_center_offset_array + (
            self.chord_length_array[0] - self.chord_length_array) / 2
        return leading_edge_x_array

    @property
    def trailing_edge_x_array(self) -> np.ndarray:
        """
        Array of x coordinates of the trailing edge.
        """
        trailing_edge_x_array = self.chord_length_array + self.x_center_offset_array + (
            self.chord_length_array[0] - self.chord_length_array) / 2
        return trailing_edge_x_array

    @property
    def center_x_array(self) -> np.ndarray:
        """
        Array of the x coordinates of the mean between leading and trailing edges coordinates.
        """
        center_x_array = (self.leading_edge_x_array +
                          self.leading_edge_x_array) / 2
        return center_x_array

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
        wing_span = np.max(self.y_array) * 2
        return wing_span

    @property
    def aspect_ratio(self):
        """
        Aspect ratio. It corresponds to the squared wing span over the reference surface.
        """
        aspect_ratio = pow(self.wing_span, 2) / self.reference_surface
        return aspect_ratio

    @property
    def sweep_angle_at_leading_edge(self):
        """
        Sweep angle at the leading edge of the wing.
        """
        sweep_angle = np.atan(
            (self.leading_edge_x_array[1] - self.leading_edge_x_array[0]) / (self.y_array[1] - self.y_array[0]))
        return sweep_angle

    @property
    def sweep_angle_at_trailing_edge(self):
        """
        Sweep angle at the trailing edge of the wing.
        """
        sweep_angle = np.atan(
            (self.trailing_edge_x_array[0] - self.trailing_edge_x_array[1]) / (self.y_array[1] - self.y_array[0]))
        return sweep_angle

    @property
    def taper_ratio(self):
        """
        Taper ratio. It corresponds to the ratio of the min chord and the max chord.
        """
        taper_ratio = np.min(self.chord_length_array) / \
            np.max(self.chord_length_array)
        return taper_ratio

    def initialize(self):
        """
        Initialize the wing by setting the undefined array to default values.

        Raises
        ------
        ValueError
            Raise error if the y array is not defined.
        """

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
        """
        Check that the wing is correctly initialized.
        """

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
        """
        Re-interpolate the wing arrays.

        Parameters
        ----------
        new_y_array : np.ndarray
            New y array on which to interpolate.
        update_y_array : bool, optional
            Indicate wether to update the y array, by default True
        """

        # Create interpolation functions
        chord_length_function = make_interp_spline(
            self.y_array, self.chord_length_array)
        twisting_angle_function = make_interp_spline(
            self.y_array, self.twisting_angle_array)
        x_center_offset_function = make_interp_spline(
            self.y_array, self.x_center_offset_array)

        # Interpolate on new points
        self.chord_length_array = chord_length_function(new_y_array)
        self.twisting_angle_array = twisting_angle_function(new_y_array)
        self.x_center_offset_array = x_center_offset_function(new_y_array)

        # Update y array if needed
        if update_y_array:
            self.y_array = new_y_array

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

        # Group to create a contour
        x_contour_array = np.concatenate(
            (self.leading_edge_x_array, self.trailing_edge_x_array[::-1], [self.leading_edge_x_array[0]]), axis=0)
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

    def create_3D_animation(self, output_path: str, nb_points_airfoil: int = 50, nb_frames: int = 60, time_step: float = 0.05, **kwargs):
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
        kwargs : dict
            Parameters to pass to the plot 3D function.
        """

        # Create the wing surface
        # wing_surface = self.create_wing_3D_surface(nb_points_airfoil)

        # pl = pv.Plotter(off_screen=True)
        # pl.add_mesh(wing_surface)

        # Prepare the pyvista scene
        pl, wing_surface = self.plot_3D(
            nb_points_airfoil=nb_points_airfoil, for_animation=True, **kwargs)

        # Animate
        pl.show(auto_close=False)
        path = pl.generate_orbital_path(
            n_points=nb_frames, shift=wing_surface.length, factor=3.0)
        if not output_path.endswith(".gif"):
            output_path = output_path + ".gif"
        pl.open_gif(output_path)
        pl.orbit_on_path(path, write_frames=True,
                         step=time_step, progress_bar=True)
        pl.close()

    def create_wing_3D_surface(self, nb_points_airfoil: int = 50):
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
        points_array = np.zeros((nb_points_airfoil * self.y_array.size, 3))
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

    def plot_3D(self,
                nb_points_airfoil: int = 50,
                show_symmetric_wing: bool = False,
                show_drag: None | Literal["blasius",
                                          "polhausen", "simulation"] = None,
                velocity_method: Literal["constant", "panels"] = "constant",
                velocity=None,
                rho=None,
                nu=None,
                for_animation=False):
        """
        Plot the shape of the wing in 3D.
        """

        # Prepare settings according to mode
        if for_animation:
            off_screen = True
        else:
            off_screen = False

        # Create the 3D surface
        wing_surface = self.create_wing_3D_surface(
            nb_points_airfoil)

        # Create the plotter object
        p = pv.Plotter(off_screen=off_screen)

        # Create a symmetry if needed
        if show_symmetric_wing:
            wing_reflected = wing_surface.reflect((0, 1, 0))
            p.add_mesh(wing_reflected)

        if show_drag is not None:
            # Raise error if not enough parameters
            if velocity is None or rho is None or nu is None:
                raise ValueError(
                    f"The values of velocity or rho or nu are not properly defined. Please provide them all to show the drag on the wing.")

            # Allocate an array for the drag
            drag_array = np.zeros(nb_points_airfoil * self.y_array.size)
            for i in range(self.y_array.size):
                points_index = i * nb_points_airfoil
                drag_array[points_index:points_index + nb_points_airfoil // 2] = self.compute_zero_lift_drag_on_wing_slice(
                    i, velocity, rho, nu, show_drag, velocity_method, nb_points_airfoil // 2, return_array=True, face="upper")[1]
                drag_array[points_index + nb_points_airfoil // 2:points_index + nb_points_airfoil] = self.compute_zero_lift_drag_on_wing_slice(
                    i, velocity, rho, nu, show_drag, velocity_method, nb_points_airfoil // 2, return_array=True, face="lower")[1][::-1]

            scalars = drag_array
            max_drag = np.max(drag_array[drag_array != np.inf])
            scalars = np.nan_to_num(scalars, nan=max_drag, posinf=max_drag)
            scalars[scalars > 1e10] = max_drag
        else:
            scalars = None

        # Plot
        p.add_mesh(wing_surface, scalars=scalars, log_scale=True)

        # Show if needed
        if not for_animation:
            p.show()
        else:
            return p, wing_surface

    def save_3D_shape(self, output_path: str, nb_points_airfoil: int = 50):
        """
        Save the wing shape as a 3D object. The format can be '.ply', '.vtp', '.stl', '.vtk', '.geo', '.obj' or '.iv'.

        Parameters
        ----------
        output_path : str
            Path of the output file.
        nb_points_airfoil : int, optional
            Number of points to use for the airfoil, by default 50
        """

        # Create the 3D surface
        wing_surface = self.create_wing_3D_surface(nb_points_airfoil)

        # Save the output file
        wing_surface.extract_geometry().save(output_path)

    def compute_fourrier_coefficients(self, alpha: float, nb_points_fourrier: int = 10):
        """
        Compute the fourrier coefficients used to determine the lift and induced drag.

        Parameters
        ----------
        alpha : float
            Angle of incidence of the wing.
        nb_points_fourrier : int, optional
            Number of points for the fourrier decomposition, by default 10

        Returns
        -------
        tuple[np.ndarray,np.ndarray]
            Tuple containing the fourrier coefficients and their ids.
        """

        # Change variable from y to theta
        theta_array = convert_y_to_theta(self.y_array[::-1], self.wing_span)

        # Create functions to be able to interpolate the arrays on new positions
        chord_length_on_theta_func = make_interp_spline(
            theta_array, self.chord_length_array[::-1])
        twisting_angle_on_theta_func = make_interp_spline(
            theta_array, self.twisting_angle_array[::-1])

        # Create theta positions to interpolate
        theta_interpolation_pos = np.linspace(
            np.pi / 2 / nb_points_fourrier, np.pi / 2, nb_points_fourrier)
        chord_length_on_theta = chord_length_on_theta_func(
            theta_interpolation_pos)
        twisting_angle_on_theta = twisting_angle_on_theta_func(
            theta_interpolation_pos)

        # Compute the alpha at zero lift for the airfoil
        alpha_zero_lift = self.base_airfoil.compute_alpha_zero_lift()

        # Allocate variables to define the system to solve
        mat_A = np.zeros((nb_points_fourrier, nb_points_fourrier))
        vec_B = np.zeros(nb_points_fourrier)

        # Iterate to fill the system
        for i in range(nb_points_fourrier):
            mu_loc = (
                2 * np.pi * chord_length_on_theta[i]) / (4 * self.wing_span)
            for j in range(nb_points_fourrier):
                n = 2 * j + 1
                mat_A[i, j] = (np.sin(theta_interpolation_pos[i]) +
                               n * mu_loc) * np.sin(n * theta_interpolation_pos[i])

            vec_B[i] = np.sin(theta_interpolation_pos[i]) * \
                (-twisting_angle_on_theta[i] +
                 alpha - alpha_zero_lift) * mu_loc

        # Solve the system
        An_vec_odd = np.linalg.solve(mat_A, vec_B)
        n_vec_odd = np.linspace(
            1, nb_points_fourrier * 2 - 1, nb_points_fourrier)

        return An_vec_odd, n_vec_odd

    def compute_lift_and_induced_drag_coefficients(self, alpha: float, nb_points_fourrier: int = 10):
        """
        Compute the coefficients of lift and induced drag.

        Parameters
        ----------
        alpha : float
            Angle of incidence of the wing.
        nb_points_fourrier : int, optional
            Number of points for the fourrier decomposition, by default 10

        Returns
        -------
        tuple[float,float]
            Tuple containing the lift and induced drag coefficients.
        """

        # Compute the fourrier coefficients
        An_vec_odd, n_vec_odd = self.compute_fourrier_coefficients(
            alpha, nb_points_fourrier)

        # Compute the lift and induced drag coefficients
        CL = np.pi * An_vec_odd[0] * self.aspect_ratio
        CD = np.pi * self.aspect_ratio * \
            np.sum(n_vec_odd * np.power(An_vec_odd, 2))

        return CL, CD

    def compute_zero_lift_drag_on_wing_slice(self,
                                             y_index: int,
                                             velocity: float,
                                             rho: float,
                                             nu: float,
                                             drag_method: Literal["blasius",
                                                                  "polhausen", "simulation"] = "polhausen",
                                             velocity_method: Literal["constant",
                                                                      "panels"] = "constant",
                                             nb_points: int = 1000,
                                             return_array: bool = False,
                                             face: Literal["upper", "lower"] = "upper"):
        """
        Compute the drag of the wing at zero lift.

        Parameters
        ----------
        y_index : int
            Index of the slice.
        velocity : float
            Velocity of the external flow.
        rho : float
            Density of the fluid.
        nu : float
            Kinematic viscosity of the fluid.
        drag_method : Literal["blasius", "polhausen", "simulation"], optional
            Method to use for the drag computation, by default "polhausen"
        velocity_method : Literal["constant", "panels"], optional
            Method to use for the velocity computation, by default "constant"
        nb_points : int, optional
            Number of points to use for the computation, by default 1000
        face : Literal["upper", "lower"]
            Indicate on which face to compute the drag.

        Returns
        -------
        float
            Drag on the wing slice.

        Raises
        ------
        NotImplementedError
            Raise error if the given method is not defined.
        """

        # Create x array
        x_array = np.linspace(
            0, self.chord_length_array[y_index], nb_points)

        # Create velocity array
        if velocity_method == "constant":
            velocity_array = velocity * np.ones(nb_points)
        elif velocity_method == "panels":
            raise NotImplementedError
            if face == "upper":
                pass
            elif face == "lower":
                pass
            else:
                raise ValueError
            velocity_array = ...
        else:
            raise NotImplementedError(
                f"The velocity method {velocity_method} is not implemented. Please use 'constant' or 'panels'.")

        # Compute the drag on the wing slice
        result = compute_linear_drag(
            x_array, velocity_array, rho, nu, drag_method, return_array=return_array)

        return result

    def compute_zero_lift_drag(self,
                               velocity: float,
                               rho: float,
                               nu: float,
                               drag_method: Literal["blasius",
                                                    "polhausen", "simulation"] = "polhausen",
                               velocity_method: Literal["constant",
                                                        "panels"] = "constant",
                               nb_points: int = 1000,
                               face: Literal["both", "upper", "lower"] = "both"):
        """
        Compute the drag of the half-wing at zero lift.

        Parameters
        ----------
        velocity : float
            Velocity of the external flow.
        rho : float
            Density of the fluid.
        nu : float
            Kinematic viscosity of the fluid.
        drag_method : Literal["blasius", "polhausen", "simulation"], optional
            Method to use for the drag computation, by default "polhausen"
        velocity_method : Literal["constant", "panels"], optional
            Method to use for the velocity computation, by default "constant"
        nb_points : int, optional
            Number of points to use for the computation, by default 1000
        face : Literal["both", "upper", "lower"]
            Indicate on which face to compute the drag.

        Returns
        -------
        float
            Drag at zero lift
        """

        if face == "both":
            zero_lift_drag = self.compute_zero_lift_drag(velocity, rho, nu, drag_method, velocity_method, nb_points, face="upper") + \
                self.compute_zero_lift_drag(
                    velocity, rho, nu, drag_method, velocity_method, nb_points, face="lower")
        else:
            # Compute linear drag on the wing
            linear_drag_on_wing = np.zeros(self.y_array.shape)
            for i in range(self.y_array.size):
                linear_drag_on_wing[i] = self.compute_zero_lift_drag_on_wing_slice(
                    i, velocity, rho, nu, drag_method, velocity_method, nb_points, face=face)

            # Integrate the drag over y
            zero_lift_drag = np.trapezoid(linear_drag_on_wing, self.y_array)

        return zero_lift_drag
