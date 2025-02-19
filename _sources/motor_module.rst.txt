Motor module
============

This module includes an electric motor model that can be used to compute rotation speed, torque efficiency based on the characteristics.

Physical background
-------------------

The model can be represented as shown below:

.. image:: ./figures/electric_motor.png
  :width: 400
  :alt: Electric diagram

It considers two resistors, one inside the motor and the other one in the power source in addition to the motor itself.

A coefficient :math:`K_v` is used to bind the rotation speed with the electromotive force :math:`E`:

.. math:: E = \frac{\omega}{K_v}
    :label: E-def

All the relations in the module derive from the electric equations behind this representation.

Database source
---------------

Some motor characteristics are based on the technical documentation of `EMRAX <https://emrax.com/e-motors/>`_. The other ones are fictional examples.
