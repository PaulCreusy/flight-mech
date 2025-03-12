"""
Module providing simple gas models.
"""

###########
# Imports #
###########

# Python imports #

from abc import ABC, abstractmethod
from pyparsing.exceptions import ParseException

# Dependencies #

import numpy as np
from periodictable.formulas import formula

# Local imports #

from flight_mech._common import lazy, hash_dict

#############
# Constants #
#############

ideal_gas_constant = 8.31446261815324
monoatomic_heat_capacity_ratio = 5 / 3
diatomic_heat_capacity_ratio = 7 / 5

# Values computed between 173K and 373K from : https://webbook.nist.gov/chemistry/fluid/#
SUTHERLAND_CONSTANTS = {
    "N2": {"T_0": 173, "mu_0": 11.415e-6, "S": 107.5},
    "H2": {"T_0": 173, "mu_0": 6.0978e-6, "S": 60.2},
    "O2": {"T_0": 173, "mu_0": 12.948e-6, "S": 124.1},
    "CH4": {"T_0": 173, "mu_0": 6.6944e-6, "S": 170.1},
    "He": {"T_0": 173, "mu_0": 13.743e-6, "S": 54.35},
}

#############
# Functions #
#############

def compute_Sutherland_constant(mu_0: float, T_0: float, mu_1: float, T_1: float):
    """
    Compute the Sutherland's constant from two known viscosities.

    Parameters
    ----------
    mu_0 : float
        Viscosity at reference temperature in kg.m-1.s-1.
    T_0 : float
        Reference temperature in K.
    mu_1 : float
        Viscosity at other temperature in kg.m-1.s-1.
    T_1 : float
        Other temperature in K.

    Returns
    -------
    float
        Sutherland's constant.
    """

    temperature_ratio = T_1 / T_0
    S = (mu_0 * np.power(temperature_ratio, 3 / 2) * T_0 - mu_1 * T_1) / \
        (mu_1 - mu_0 * np.power(temperature_ratio, 3 / 2))

    return S

def compute_viscosity_Sutherland(mu_0: float, T_0: float, S: float, T_1: float):
    """
    Compute the viscosity of a gas knowing the Sutherland's constant and the viscosity
    at a reference temperature.

    Parameters
    ----------
    mu_0 : float
        Viscosity at reference temperature in kg.m-1.s-1.
    T_0 : float
        Reference temperature in K.
    S : float
        Sutherland's constant.
    T_1 : float
        Temperature in K.

    Returns
    -------
    float
        Viscosity at target temperature in kg.m-1.s-1.
    """

    temperature_ratio = T_1 / T_0
    mu = mu_0 * np.power(temperature_ratio, 3 / 2) * \
        (T_0 + S) / (T_1 + S)

    return mu

def compute_mixing_viscosity_Wilke(molar_fraction_array: np.ndarray, viscosity_array: np.ndarray, molar_mass_array: np.ndarray):
    """
    Compute the viscosity of a gas mixture using Wilke's equation.

    Reference
    ---------
    Based on: A viscosity for Gas Mixture, C. R. Wilke, 1950, The Journal of Chemical Physics, Volume 18, Number 4

    Parameters
    ----------
    molar_fraction_array : np.ndarray
        Array containing the molar fraction of each gas.
    viscosity_array : np.ndarray
        Array containing the viscosity of each gas.
    molar_mass_array : np.ndarray
        Array containing the molar mass of each gas.

    Returns
    -------
    float
        Viscosity of the mixture.
    """

    # Extract the number of gases
    nb_gases = molar_fraction_array.shape[0]

    # Compute the phi matrix
    phi_matrix = np.power(1 + np.power(viscosity_array.reshape(nb_gases, 1) / viscosity_array.reshape(
        1, nb_gases), 0.5) * np.power(molar_mass_array.reshape(1, nb_gases) / molar_mass_array.reshape(nb_gases, 1), 0.5), 2) / \
        ((4 / np.sqrt(2)) *
         np.power(1 + (molar_mass_array.reshape(nb_gases, 1) /
                       molar_mass_array.reshape(1, nb_gases)), 0.5))

    # Compute the viscosity mixture
    mixture_viscosity = np.sum(viscosity_array / (np.dot(phi_matrix,
                               molar_fraction_array.reshape(nb_gases, 1)).flatten() / molar_fraction_array))

    return mixture_viscosity

###########
# Classes #
###########

class GasModel(ABC):
    """
    Gas abstract class.
    """

    name: str
    temperature: float = 273  # K
    pressure: float = 1e5  # Pa
    molar_mass: float = 1.  # kg.mol-1

    def __init__(self, name: str, molar_mass: float | None = None, temperature: float | None = None, pressure: float | None = None):
        self.name = name

        # Attribute molar mass
        if molar_mass is None:
            # Check if name can be interpreted as a valid formula
            try:
                self.molar_mass = formula(name).mass * 1e-3
            except ParseException as e:
                raise Warning(
                    f"The given name ({name}) is not a valid chemical formula. A molar mass must be provided.")
        else:
            self.molar_mass = molar_mass

        if temperature is not None:
            self.temperature = temperature

        if pressure is not None:
            self.pressure = pressure

    def __hash__(self):
        data_dict = {}
        for key in dir(self):
            value = getattr(self, key)
            if isinstance(value, float):
                data_dict[key] = value
        return hash_dict(data_dict)

    @property
    @abstractmethod
    def density(self) -> float:
        """
        Density of the gas in kg.m-3.
        """
        pass

    @property
    @abstractmethod
    def Cp(self) -> float:
        """
        Heat capacity at constant pressure in J.K-1.
        """
        pass

    @property
    @abstractmethod
    def Cv(self) -> float:
        """
        Heat capacity at constant volume in J.K-1.
        """
        pass

    @property
    @abstractmethod
    def sound_velocity(self) -> float:
        """
        Sound velocity in m.s-1.
        """
        pass

    @property
    @abstractmethod
    def dynamic_viscosity(self) -> float:
        """
        Dynamic viscosity in kg.m-1.s-1.
        """
        pass

    @property
    def kinematic_viscosity(self):
        """
        Dynamic viscosity in m2.s-1.
        """
        return self.dynamic_viscosity / self.rho

    @property
    def rho(self):
        """
        Alias for density.
        """
        return self.density

    @property
    def mu(self):
        """
        Alias for dynamic viscosity.
        """
        return self.dynamic_viscosity

    @property
    def nu(self):
        """
        Alias for kinematic viscosity.
        """
        return self.kinematic_viscosity

    @property
    def P(self) -> float:
        """
        Alias for pressure.
        """
        return self.pressure

    @P.setter
    def P(self, value: float):
        self.pressure = value

    @property
    def T(self) -> float:
        """
        Alias for temperature.
        """
        return self.temperature

    @T.setter
    def T(self, value: float):
        self.temperature = value


class PerfectGas(GasModel):
    """
    Class to define a perfect gas and its properties.
    """

    heat_capacity_ratio: float
    Sutherland_constant: float | None = None
    mu_0: float | None = None
    T_0: float | None = None

    @property
    def gamma(self):
        """
        Alias for heat capacity ratio.
        """
        return self.heat_capacity_ratio

    @gamma.setter
    def gamma(self, value: float):
        self.heat_capacity_ratio = value

    @property
    def density(self):
        rho = self.pressure / (self.r * self.temperature)
        return rho

    @property
    def Cp(self):
        Cp = self.gamma / (self.gamma - 1) * self.r
        return Cp

    @property
    def Cv(self):
        Cv = (1 / (self.gamma - 1)) * self.r
        return Cv

    @property
    def sound_velocity(self):
        sound_velocity = np.sqrt(self.gamma * self.r * self.temperature)
        return sound_velocity

    @property
    def dynamic_viscosity(self):
        if self.Sutherland_constant is not None and self.mu_0 is not None and self.T_0 is not None:
            temperature_ratio = self.T / self.T_0
            mu = self.mu_0 * np.power(temperature_ratio, 3 / 2) * (
                self.T_0 + self.Sutherland_constant) / (self.T + self.Sutherland_constant)
            return mu

        if self.name in SUTHERLAND_CONSTANTS:
            self.Sutherland_constant = SUTHERLAND_CONSTANTS[self.name]["S"]
            self.mu_0 = SUTHERLAND_CONSTANTS[self.name]["mu_0"]
            self.T_0 = SUTHERLAND_CONSTANTS[self.name]["T_0"]
            return self.dynamic_viscosity

        raise ValueError(
            "No values available to compute the dynamic viscosity. Please provide a Sutherland's constant and a reference viscosity and temperature.")

    @property
    def specific_gas_constant(self):
        """
        Specific gas constant in J.kg-1.K-1.
        """
        specific_gas_constant = ideal_gas_constant / self.molar_mass
        return specific_gas_constant

    @property
    def r(self):
        """
        Alias for specific gas constant.
        """
        return self.specific_gas_constant


class MonoatomicPerfectGas(PerfectGas):
    """
    Class to define a monoatomic perfect gas and its properties.
    """

    heat_capacity_ratio = monoatomic_heat_capacity_ratio

class DiatomicPerfectGas(PerfectGas):
    """
    Class to define a diatomic perfect gas and its properties.
    """

    heat_capacity_ratio = diatomic_heat_capacity_ratio

class GasMixture(GasModel):
    """
    Class to define a gas mixture.
    """

    gas_model_dict: dict[str, GasModel]
    gas_molar_fraction_dict: dict[str, float]
    _gas_state_hash: int = 0
    _molar_mass: float

    def __init__(self,
                 gas_model_dict: dict[str, GasModel],
                 gas_molar_fraction_dict: dict[str, float],
                 temperature: float | None = None,
                 pressure: float | None = None):

        self.gas_model_dict = gas_model_dict
        self.gas_molar_fraction_dict = gas_molar_fraction_dict
        self._check_gaz_molar_fraction()
        self._check_dict_correspondance()

        if temperature is not None:
            self.temperature = temperature

        if pressure is not None:
            self.pressure = pressure

    def __hash__(self):
        return self._get_gas_state_hash() + hash(self.T) + hash(self.P)

    def _check_gaz_molar_fraction(self):
        """
        Verify that the gas molar fraction dict is correctly defined.

        Raises
        ------
        ValueError
            Raise error if sum is not equal to 1 or in case of negative values or greater than one.
        """

        total_molar_fraction = 0.
        for key in self.gas_molar_fraction_dict:
            value = self.gas_molar_fraction_dict[key]
            if not (0 <= value <= 1):
                raise ValueError(
                    f"The fraction of the gas {key} is outside interval [0;1] ({value}).")
            total_molar_fraction += value
        if not np.isclose(total_molar_fraction, 1.):
            raise ValueError(
                f"The total fraction is not equal to 1. ({total_molar_fraction}).")

    def _check_dict_correspondance(self):
        """
        Verify that the model and fraction dict contain the same keys.

        Raises
        ------
        ValueError
            Raise error if the keys are different.
        """

        if sorted(list(self.gas_model_dict)) != sorted(
                list(self.gas_molar_fraction_dict)):
            raise ValueError(
                "The keys of the model and fraction dictionaries are not the same.")

    def _update_pressure(self):
        """
        Update the pressure variable of all gas models
        """
        for key in self.gas_model_dict:
            self.gas_model_dict[key].pressure = self.pressure

    def _update_temperature(self):
        """
        Update the temperature variable of all gas models
        """
        for key in self.gas_model_dict:
            self.gas_model_dict[key].temperature = self.temperature

    def _update_gas_variables(self):
        """
        Update the variables of all gas models.
        """
        self._update_pressure()
        self._update_temperature()

    def _get_gas_state_hash(self):
        """
        Get the hash of the gas state.
        """
        gas_state_hash = hash_dict(self.gas_model_dict) +\
            hash_dict(self.gas_molar_fraction_dict)
        return gas_state_hash

    def _update_gas_state_hash(self):
        """
        Update the hash of the gas state.
        """
        self._gas_state_hash = self._get_gas_state_hash()

    @property
    def _is_gas_model_dict_hash_correct(self):
        return self._gas_state_hash == self._get_gas_state_hash()

    @property
    @lazy
    def molar_mass(self) -> float:
        """
        Molar mass of the mixture in kg.mol-1.
        """

        # Update gas variables
        self._update_gas_variables()

        # Compute molar mass
        molar_mass = 0
        for key in self.gas_model_dict:
            molar_mass += self.gas_molar_fraction_dict[key] * \
                self.gas_model_dict[key].molar_mass

        return molar_mass

    @property
    @lazy
    def gas_mass_fraction_dict(self) -> dict[str, float]:
        """
        Dictionary containing the mass fraction of each gas in the mixture.
        """

        # Update gas variables
        self._update_gas_variables()

        gas_mass_fraction_dict = {}
        for key in self.gas_model_dict:
            gas_mass_fraction_dict[key] = self.gas_model_dict[key].molar_mass * self.gas_molar_fraction_dict[key] / \
                self.molar_mass

        return gas_mass_fraction_dict

    @property
    @lazy
    def density(self):

        # Update gas variables
        self._update_gas_variables()

        density_invert = 0
        for key in self.gas_model_dict:
            density_invert += self.gas_mass_fraction_dict[key] / \
                self.gas_model_dict[key].density
        density = 1 / density_invert

        return density

    @property
    @lazy
    def Cp(self):

        # Update gas variables
        self._update_gas_variables()

        Cp = 0
        for key in self.gas_model_dict:
            Cp += self.gas_model_dict[key].Cp * \
                self.gas_molar_fraction_dict[key]
        return Cp

    @property
    @lazy
    def Cv(self):

        # Update gas variables
        self._update_gas_variables()

        Cv = 0
        for key in self.gas_model_dict:
            Cv += self.gas_model_dict[key].Cv * \
                self.gas_molar_fraction_dict[key]
        return Cv

    @property
    @lazy
    def sound_velocity(self):

        # Update gas variables
        self._update_gas_variables()

        sound_velocity_invert = 0
        for key in self.gas_model_dict:
            sound_velocity_invert += self.gas_molar_fraction_dict[key] / \
                self.gas_model_dict[key].sound_velocity
        sound_velocity = 1 / sound_velocity_invert
        return sound_velocity

    @property
    @lazy
    def dynamic_viscosity(self):

        # Update gas variables
        self._update_gas_variables()

        # Extract arrays of values
        molar_fraction_array = np.array(
            [self.gas_molar_fraction_dict[key] for key in self.gas_model_dict])
        molar_mass_array = np.array(
            [self.gas_model_dict[key].molar_mass for key in self.gas_model_dict])
        viscosity_array = np.array(
            [self.gas_model_dict[key].dynamic_viscosity for key in self.gas_model_dict])

        # Compute mixture viscosity
        mixture_viscosity = compute_mixing_viscosity_Wilke(
            molar_fraction_array, viscosity_array, molar_mass_array)

        return mixture_viscosity

class Air(GasMixture):

    def __init__(self, temperature: float | None = None, pressure: float | None = None):

        gas_model_dict = {
            "O2": DiatomicPerfectGas("O2"),
            "N2": DiatomicPerfectGas("N2")
        }
        gas_molar_fraction_dict = {
            "O2": 0.21,
            "N2": 0.79
        }

        super().__init__(gas_model_dict, gas_molar_fraction_dict, temperature, pressure)
