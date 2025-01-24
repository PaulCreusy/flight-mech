Atmosphere module
=================

This module defines a simple atmospheric model that is mostly used to compute air density. The atmospheric model used is a simplification of the standard atmosphere and it is valid between 0 and 11km of altitude.

Physical background
-------------------

The model is based here on a single formula:

.. math:: \sigma = \frac{20-z}{20+z}
   :label: sigma-formula

with:

.. math:: \sigma = \frac{\rho}{\rho_0}
    :label: sigma-def

Examples
--------

The functions allow to compute the density:

>>> compute_sigma_from_altitude(4.2e3)
0.6528925619834711

or directly the mass per cubic meter:

>>> compute_air_density_from_altitude(3e3)
0.9054347826086957

For a given sigma, it is also possible to compute the altitude:

>>> compute_altitude_from_sigma(0.7891)
2357.609971494047