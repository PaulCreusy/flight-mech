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
    def rho(self):
        """
        Alias for density.
        """
        return self.density

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

    def __init__(self, gas_model_dict: dict[str, GasModel], gas_molar_fraction_dict: dict[str, float]):
        self.gas_model_dict = gas_model_dict
        self.gas_molar_fraction_dict = gas_molar_fraction_dict
        self._check_gaz_molar_fraction()
        self._check_dict_correspondance()

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
