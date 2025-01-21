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
    x_array: np.ndarray

    def __init__(self):
        pass

    @property
    def camber_y_array(self) -> np.ndarray:
        camber_y_array = (self.extrados_y_array + self.intrados_y_array) / 2
        return camber_y_array

    @property
    def max_thickness(self):
        pass

    @property
    def max_camber(self):
        pass

    def load_database_airfoil(self, airfoil_name: str, airfoil_data_folder: str = default_airfoil_database):
        file_path = os.path.join(airfoil_data_folder, airfoil_name + ".txt")
        self.load_selig_file(file_path, skiprows=1)

    def load_selig_file(self, file_path: str, skiprows: int = 0):
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

    def plot(self):
        x_plot = np.concatenate(
            (self.x_array, self.x_array[::-1], [self.x_array[0]]), axis=0)
        y_plot = np.concatenate(
            (self.extrados_y_array, self.intrados_y_array[::-1], [self.extrados_y_array[0]]), axis=0)
        plt.plot(x_plot, y_plot)
        plt.axis("equal")
        plt.show()

    def re_interpolate(new_x_array):
        pass

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

        try:
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
        except requests.exceptions.RequestException as e:
            raise e
