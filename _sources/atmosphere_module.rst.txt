Atmosphere module
=================

This module defines several atmosphere models to compute density, pressure, temperature and other useful quantities at the altitude of interest for the other modules. Here is the list of the atmosphere models implemented:

* Constant atmosphere
* Linear atmosphere
* Standard atmosphere

Constant atmosphere
-------------------

This is a debug atmosphere model. The density is considered constant equal to the sea level at any altitude.

Linear atmosphere
-----------------

This is a simple atmospheric model that is mostly used to compute air density. The atmospheric model used is a simplification of the standard atmosphere and it is valid between 0 and 11km of altitude.

Physical background
^^^^^^^^^^^^^^^^^^^

The model is based here on a single formula:

.. math:: \sigma = \frac{20-z}{20+z}
   :label: sigma-formula

with:

.. math:: \sigma = \frac{\rho}{\rho_0}
    :label: sigma-def

Examples
^^^^^^^^

The functions allow to compute the density:

>>> compute_sigma_from_altitude(4.2e3)
0.6528925619834711

or directly the mass per cubic meter:

>>> compute_density_from_altitude(3e3)
0.9054347826086957

For a given sigma, it is also possible to compute the altitude:

>>> compute_altitude_from_sigma(0.7891)
2357.609971494047

Standard atmosphere
-------------------

This is the International Standard Atmosphere model and it allows to compute the following quantities:

* temperature
* pressure
* density
* sound speed
* dynamic viscosity
* kinematic viscosity

Physical background
^^^^^^^^^^^^^^^^^^^

Here is the list of hypotheses of the model:

1. The air is considered to be a perfect gas with :math:`r = 287 J.kg^{-1}.K`.
2. The air is dry.
3. The gravity force is uniform.
4. The atmosphere is at equilibrium. Therefore :math:`dp = -\rho g_0 dh`
5. The temperature changes with altitude following these relations:

    1. Between 0 and 11km: :math:`T(h) = 288.15 + h (216.65 - 288.15)/11000`
    2. Between 11 and 20km: :math:`T(h) = 216.65`
    3. Between 20 and 32km: :math:`T(h) = 216.65 + (h - 20000) (228.65 - 216.65)/12000`
    4. Between 32 and 47km: :math:`T(h) = 228.65 + (h - 32000) (270.65 - 228.65)/12000`

The formulas of the other quantities directly derive from theses hypotheses.
