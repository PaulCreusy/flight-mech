Gas Module
==========

This module provides simple gas models for thermodynamic calculations. The models implemented include perfect gases and gas mixtures, each with specific assumptions and limitations.

.. module:: gas_models
   :synopsis: Provides simple gas models for other modules.


GasModel (Abstract Base Class)
------------------------------

The `GasModel` serves as an abstract base class for defining different gas models. It provides a common structure for all gases. It takes two quantities in input, the temperature (:math:`T`) and the pressure (:math:`P`) and uses them to compute the following properties:

* :math:`\rho` : the density in :math:`kg.m^{-3}`
* :math:`C_p` : the heat capacity at constant pressure in :math:`J.K^{-1}`
* :math:`C_v` : the heat capacity at constant volume in :math:`J.K^{-1}`
* :math:`c` : the sound velocity in :math:`m.s^{-1}`

**Assumptions:**

* Represents a generic gas with thermodynamic properties.
* Assumes thermodynamic equilibrium.

**Limitations:**

* Requires subclasses to define actual gas properties.

.. hint:: 
   If you pass a valid chemical formula (such as "O2") as name when defining a gas, its molar mass will be computed automatically.

Perfect Gas Model
-----------------

The `PerfectGas` model follows the ideal gas law:

.. math:: PV = nRT

where:

* :math:`P` : the pressure,
* :math:`V` : the volume,
* :math:`n` : the number of moles,
* :math:`R` : the universal gas constant,
* :math:`T` : the temperature.

The other quantities are defined as follows:

.. math:: C_p = \frac{\gamma}{\gamma - 1} r

.. math:: C_v = \frac{1}{\gamma - 1} r

.. math:: c = \sqrt{\gamma r T}

**Assumptions:**

* The gas behaves ideally with no intermolecular forces.
* The heat capacity ratio (:math:`\gamma`) remains constant.

**Limitations:**

* Does not account for real gas effects such as compressibility at high pressures.
* Assumes uniform temperature and pressure.

For a `MonoatomicPerfectGas`, the heat capacity ratio is fixed at 5/3.

For a `DiatomicPerfectGas`, the heat capacity ratio is set to 7/5.

Gas Mixture Model
-----------------

The `GasMixture` model represents a mixture of gases, where the overall properties are derived from the weighted contributions of individual gas components. It assumes Dalton’s Law of Partial Pressures and perfect mixing:

.. math:: P_{total} = \sum P_i

where :math:`P_i` is the partial pressure of each gas.

The molar mass of the mixture is computed as:

.. math:: M_{mix} = \sum x_i M_i

where:

* :math:`x_i` is the molar fraction of each gas,
* :math:`M_i` is the molar mass of each gas.

The heat capacity of the mixture is computed as a weighted sum of each gas:

.. math:: Cp_{mix} = \sum x_i Cp_i

.. math:: Cv_{mix} = \sum x_i Cv_i

The sound speed velocity and density are computed using the inverse molar fraction weighting rule:

.. math:: \frac{1}{\rho_{mix}} = \sum \frac{x_i}{\rho_i}

.. math:: \frac{1}{c_{mix}} = \sum \frac{x_i}{c_i}

Other properties, such as density and heat capacities, are obtained as weighted sums of the individual gas contributions.

**Assumptions:**

* The gases mix ideally without chemical interactions.
* Follows Dalton’s Law for partial pressures.
* Properties are computed as weighted averages.

**Limitations:**

* Neglects non-ideal behavior and interactions between different gas species.

Usage Example
-------------

Here is a short example on how to use the module. Please check the :ref:`/examples.rst` section for more details.

.. code-block:: python

   from gas_models import MonoatomicPerfectGas

   helium = MonoatomicPerfectGas("He")
   helium.pressure = 1e5  # Pa
   helium.temperature = 273  # K

   print(helium.density)
   print(helium.Cp)
