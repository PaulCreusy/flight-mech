Flight-Mech Documentation
=========================

.. meta::
    :description lang=en:
        A Python library to design and evaluate the performances of a plane based on flight mechanics equations.

Flight-Mech includes a set of modules based on simple physical models to design a plane or compute its characteristics.

Introduction
------------

The goal of flight-mech is to describe with simple physical models, all the main components that can be defined to design a plane and ultimately with a few key dimensions and properties to evaluate the behavior of the plane.

Flight-Mech is for you if:

* You want to design a plane and are looking for a way to determine its main characteristics.
* You do not have a precise geometry or mesh.
* You want to optimize a plane over a large set of parameters.

Flight-Mech is probably not the best solution if:

* You need to evaluate the aerodynamics of a specific geometry.
* You want to simulate precisely the behavior of a plane or a component.
* Your plane can fly at supersonic speed or uses cutting-edge technologies.

Organisation
------------

This package is decomposed in various independent modules, each one describing a physical aspect of the plane. The modules implemented are the following:

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

You can install flight-mech using Pip with:

.. code:: console

   pip install flight-mech

The source code of flight-mech is available on Github at this `link <https://github.com/PaulCreusy/flight-mech/tree/main>`_.


.. toctree::
   :maxdepth: 1
   :caption: User guide

   atmosphere_module
   plane_module
   airfoil_module
   wing_module
   turbine_module
   motor_module

.. toctree:: 
   :maxdepth: 1
   :caption: Examples

.. toctree::
   :maxdepth: 3
   :caption: API Reference

   flight_mech
