Flight-mech documentation
=========================

This library includes a simple atmospheric model and a set of flight mechanics equations allowing to compute plane and airfoils characteristics.

Introduction
------------

This software includes various modules to build a numerical plane model and compute its characteristics. The modules implemented are the following:

* **atmosphere** : defines several atmosphere models to compute density, temperature, pressure and other quantities.
* **aerodynamics** : contains functions to compute quantities in the boundary layer of a fluid flow.
* **airfoil** : allows to define the geometry of an airfoil and compute the lift and moment coefficients.
* **wing** : allows to define the geometry of a wing and compute the lift and drag coefficients.
* **fuel** : defines several types of broadly used fuels in aeronautics.
* **turbine** : allows to define several types of turbine to compute their thrust and consumption at various operating conditions.
* **motor** : allows to define an electric motor model and compute its outputs.
* **plane** : allows to define a numerical plane model, binding the previous modules, to compute its flight characteristics. 

.. Warning::
   All units are defined in the international unit system except if explicitly indicated.

Installation
------------

You can install it using Pip with:

.. code-block:: bash

   pip install flight-mech

The source code is available on Github at this `link <https://github.com/PaulCreusy/flight-mech/tree/main>`_.


.. toctree::
   :maxdepth: 1
   :caption: Getting started

   atmosphere_module
   plane_module
   airfoil_module
   wing_module
   turbine_module
   motor_module

.. toctree::
   :maxdepth: 3
   :caption: Package documentation

   flight_mech
