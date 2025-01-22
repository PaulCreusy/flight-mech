"""
Module to analyse airfoil aerodynamic properties.
"""

###########
# Imports #
###########

# Python imports #

import os
from typing import Literal
import json

# Dependencies #

import numpy as np
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt
import pandas as pd
import requests
from bs4 import BeautifulSoup

#############
# Constants #
#############

# Define a default plane database location
default_airfoil_database = os.path.join(
    os.path.dirname(__file__), "airfoil_database")

#############
# Functions #
#############

def convert_float_or_none_to_string(value: float | None):
    if value is None:
        return ""
    return str(value)

def download_web_file(url, output_file):
    # Make the GET request to download the file
    response = requests.get(url)
    response.raise_for_status()

    # Write the file content to the specified output file
    with open(output_file, 'wb') as file:
        file.write(response.content)


###########
# Classes #
###########

class Airfoil:

    extrados_y_array: np.ndarray
    intrados_y_array: np.ndarray
    name: str
    _x_array: np.ndarray
    _chord_length: float

    def __init__(self, airfoil: str | None = None):
        if airfoil in self.list_airfoils_in_database():
            self.load_database_airfoil(airfoil)
        elif airfoil is not None:
            self.load_selig_file(airfoil)

    @property
    def camber_y_array(self) -> np.ndarray:
        camber_y_array = (self.extrados_y_array + self.intrados_y_array) / 2
        return camber_y_array

    @property
    def max_thickness(self) -> float:
        max_thickness = np.max(self.extrados_y_array -
                               self.intrados_y_array) / self.chord_length
        return max_thickness

    @property
    def max_thickness_location(self) -> float:
        max_thickness_idx = np.argmax(self.extrados_y_array -
                                      self.intrados_y_array)
        max_thickness_location = self.x_array[max_thickness_idx]
        return max_thickness_location

    @property
    def max_camber(self) -> float:
        max_camber = np.max(np.abs(self.camber_y_array -
                            self.chord_y_array)) / self.chord_length
        return max_camber

    @property
    def max_camber_location(self) -> float:
        max_camber_idx = np.argmax(np.abs(self.camber_y_array -
                                          self.chord_y_array))
        max_camber_location = self.x_array[max_camber_idx]
        return max_camber_location

    @property
    def x_array(self) -> np.ndarray:
        return self._x_array

    @x_array.setter
    def x_array(self, value: np.ndarray):
        self._x_array = value
        self._chord_length = np.max(value)

    @property
    def chord_length(self) -> float:
        return self._chord_length

    @chord_length.setter
    def chord_length(self, value):
        ratio = value / self._chord_length
        self.x_array = self.x_array * ratio
        self.extrados_y_array = self.extrados_y_array * ratio
        self.intrados_y_array = self.intrados_y_array * ratio

    @property
    def chord_y_array(self):
        x_0 = 0
        x_1 = self.chord_length
        y_0 = self.extrados_y_array[0]
        y_1 = self.extrados_y_array[-1]
        slope = (y_1 - y_0) / (x_1 - x_0)
        return slope * self.x_array + y_0

    def list_airfoils_in_database(self, airfoil_data_folder: str = default_airfoil_database):
        """
        Return the list of airfoils stored in the database.

        Parameters
        ----------
        airfoil_data_folder : str, optional
            Name of the airfoil database folder, by default default_airfoil_database

        Returns
        -------
        list
            List of airfoils stored in the database.
        """

        file_names_list = os.listdir(airfoil_data_folder)
        airfoil_names_list = [e.replace(".txt", "") for e in file_names_list]

        return airfoil_names_list

    def load_database_airfoil(self, airfoil_name: str, airfoil_data_folder: str = default_airfoil_database):
        """
        Load an airfoil contained in the database.

        Parameters
        ----------
        airfoil_name : str
            Name of the airfoil.
        airfoil_data_folder : str, optional
            Folder containing the airfoil file, by default default_airfoil_database
        """

        file_path = os.path.join(airfoil_data_folder, airfoil_name + ".txt")
        self.load_selig_file(file_path, skiprows=1)

    def load_selig_file(self, file_path: str, skiprows: int = 1):
        """
        Load an airfoil contained in a selig txt file.

        Parameters
        ----------
        file_path : str
            Path of the file to load.
        skiprows : int, optional
            Number of rows to skip at the beginning of the file, by default 1
        """

        # Extract airfoil name
        with open(file_path, "r") as file:
            first_line = file.readline()
        self.name = first_line.replace("\n", "")

        # Extract airfoil coordinates
        file_content = np.loadtxt(file_path, skiprows=skiprows)
        x_selig_array = file_content[:, 0]
        y_selig_array = file_content[:, 1]
        self.import_xy_selig_arrays(x_selig_array, y_selig_array)

    def import_xy_selig_arrays(self, x_selig_array: np.ndarray, y_selig_array: np.ndarray):
        # Compute the x derivative to split parts of the airfoil
        x_diff = np.zeros(x_selig_array.shape)
        x_diff[:-1] = x_selig_array[1:] - x_selig_array[:-1]
        x_diff[-1] = x_diff[-2]

        # Split parts
        diff_sign = x_diff * x_diff[0]
        common_point = np.argmax(diff_sign < 0)
        diff_sign[common_point] = 0
        part_1_x = x_selig_array[diff_sign >= 0]
        part_1_y = y_selig_array[diff_sign >= 0]
        part_2_x = x_selig_array[diff_sign <= 0]
        part_2_y = y_selig_array[diff_sign <= 0]

        # Extract all x locations
        self.x_array = np.unique(x_selig_array)

        # Reorder arrays for interpolation
        part_1_order = np.argsort(part_1_x)
        part_2_order = np.argsort(part_2_x)

        # Interpolate both parts on all x locations
        part_1_y_interpolated = np.interp(
            self.x_array, part_1_x[part_1_order], part_1_y[part_1_order])
        part_2_y_interpolated = np.interp(
            self.x_array, part_2_x[part_2_order], part_2_y[part_2_order])

        # Assign depending on which part is extrados or intrados
        if np.mean(part_1_y_interpolated) > np.mean(part_2_y_interpolated):
            self.extrados_y_array = part_1_y_interpolated
            self.intrados_y_array = part_2_y_interpolated
        else:
            self.extrados_y_array = part_2_y_interpolated
            self.intrados_y_array = part_1_y_interpolated

    def plot(self, hold_plot=False, show_chord=False, show_camber_line=False):
        """
        Plot the geometry of the airfoil.
        """

        # Concatenate extrados and intrados to plot a single line
        x_plot = np.concatenate(
            (self.x_array, self.x_array[::-1], [self.x_array[0]]), axis=0)
        y_plot = np.concatenate(
            (self.extrados_y_array, self.intrados_y_array[::-1], [self.extrados_y_array[0]]), axis=0)

        # Plot extrados and intrados
        plt.plot(x_plot, y_plot, label=self.name)
        plt.axis("equal")

        # Plot chord if needed
        if show_chord:
            plt.plot(self.x_array, self.chord_y_array,
                     "--", label=f"chord {self.name}")

        # Plot camber if needed
        if show_camber_line:
            plt.plot(self.x_array, self.camber_y_array,
                     "--", label=f"camber {self.name}")

        # Display the plot if needed
        if hold_plot is False:
            plt.legend()
            plt.title("Airfoil visualization")
            plt.show()

    def re_interpolate(self, new_x_array: np.ndarray):
        """
        Re-interpolate the extrados and intrados on the given array.

        Parameters
        ----------
        new_x_array : np.ndarray
            New x array on which to interpolate.
        """

        extrados_function = make_interp_spline(
            self.x_array, self.extrados_y_array)
        intrados_function = make_interp_spline(
            self.x_array, self.intrados_y_array)
        self.extrados_y_array = extrados_function(new_x_array)
        self.intrados_y_array = intrados_function(new_x_array)
        self.x_array = new_x_array
        self.chord_length = ...

    def import_from_airfoiltools(
            self,
            airfoil_name: str = "",
            max_thickness: float | None = None,
            min_thickness: float | None = None,
            max_camber: float | None = None,
            min_camber: float | None = None,
            maximise_glide_ratio_at_reynolds: Literal["50k", "100k",
                                                      "200k", "500k", "1M", "2M", "5M"] | None = None,
            airfoil_data_folder: str = default_airfoil_database):
        """
        Import an airfoil from airfoiltools.

        Parameters
        ----------
        airfoil_name : str, optional
            Name of the airfoil to search, by default ""
        max_thickness : float | None, optional
            Max thickness value in percent, by default None
        min_thickness : float | None, optional
            Min thickness value in percent, by default None
        max_camber : float | None, optional
            Max camber value in percent, by default None
        min_camber : float | None, optional
            Min camber value in percent, by default None
        maximise_glide_ratio_at_reynolds : Literal[&quot;50k&quot;, &quot;100k&quot;, &quot;200k&quot;, &quot;500k&quot;, &quot;1M&quot;, &quot;2M&quot;, &quot;5M&quot;] | None, optional
            Indicate Reynolds number to sort by optimal ratio, by default None
        airfoil_data_folder : str, optional
            Folder containing the airfoil database, by default default_airfoil_database

        Warning
        -------
        The thickness and camber values must be expressed in percents.

        Raises
        ------
        ValueError
            Raise error if no corresponding profile is found.
        """

        # Define the URL for the search
        search_url = "http://airfoiltools.com/search/index"

        # Set the sort mode
        if maximise_glide_ratio_at_reynolds is not None:
            sort_mode = str(9 + ["50k", "100k", "200k", "500k", "1M",
                                 "2M", "5M"].index(maximise_glide_ratio_at_reynolds))
        else:
            sort_mode = "1"

        # Convert input parameters to string
        max_thickness = convert_float_or_none_to_string(max_thickness)
        min_thickness = convert_float_or_none_to_string(min_thickness)
        max_camber = convert_float_or_none_to_string(max_camber)
        min_camber = convert_float_or_none_to_string(min_camber)

        # Set the query parameters
        params = {
            "MAirfoilSearchForm[textSearch]": airfoil_name,
            "MAirfoilSearchForm[maxThickness]": max_thickness,
            "MAirfoilSearchForm[minThickness]": min_thickness,
            "MAirfoilSearchForm[maxCamber]": max_camber,
            "MAirfoilSearchForm[minCamber]": min_camber,
            "MAirfoilSearchForm[grp]": "",
            "MAirfoilSearchForm[sort]": sort_mode,
            "yt0": "Search"
        }

        # Make the GET request to the search page
        response = requests.get(search_url, params=params)
        response.raise_for_status()

        # Parse the HTML content
        soup = BeautifulSoup(response.text, 'html.parser')

        # Find the first research proposal link (specific to the provided element structure)
        proposal_link = soup.find('a', string="Airfoil details")

        if proposal_link and 'href' in proposal_link.attrs:
            airfoil_code = proposal_link["href"].split("?")[1]
            selig_file_link = "http://airfoiltools.com/airfoil/seligdatfile?" + airfoil_code
            airfoil_name = airfoil_code.replace(
                "airfoil=", "").replace("-il", "")
            output_file_path = os.path.join(
                airfoil_data_folder, airfoil_name + ".txt")
            download_web_file(selig_file_link, output_file_path)
            print(f"Airfoil {airfoil_name} successfully downloaded.")
            self.load_database_airfoil(airfoil_name, airfoil_data_folder)
        else:
            raise ValueError("No corresponding airfoil found.")
