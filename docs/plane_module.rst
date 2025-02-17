Plane module
============

This module includes a set of flight mechanics equations allowing to compute the plane characteristics.

**Please note that all equations and variables are defined in the international unit system.**

Functionalities
---------------

The plane model allows to compute the following quantities:

* max glide ratio
* speed at specific angle of incidence and altitude
* drag
* lift
* thrust
* stall speed
* reference speed
* minimum descent gliding slope
* gliding speed
* maximum gliding time
* maximum gliding range
* authorized velocity interval at fixed thrust for flight at constant altitude
* thrust needed at fixed altitude and angle of incidence
* minimum thrust needed at fixed altitude
* speed at minimum thrust
* maximum flight altitude
* speed for maximum ascension speed
* ascension slope for a specific angle of incidence and altitude
* load factor in turn
* maximum range at fixed altitude
* maximum range at fixed speed
* endurance
* take off distance without friction
* take off distance with friction
* landing distance
* take off speed
* landing speed
* alpha and delta coefficient at a flight point

Additionally, the following graphs can be generated:

* polar graph
* thrust-speed graph
* power-speed graph

Some examples are provided in the `examples` folder (please note that they do not cover all the use cases) as well with a few plane models in the `plane_database` folder. 

Quick start
-----------

When loading a plane from the database, you can directly compute its basic characteristics like the glide ratio for instance:

>>> plane = Plane("cessna_172", "./plane_database")
>>> print(plane.C_L_f_max, plane.f_max)
0.7745966692414834, 12.909944487358056

If some information are missing, you can set them yourself.

>>> plane.m_fuel = 136.26  # kg
>>> plane.update_variables(force=True)

The module also includes functions to compute quantities at a specific flight point:

>>> plane.compute_reference_speed(8000) # m.s-1
56.214394963985406

>>> plane.compute_velocity_interval_for_fixed_thrust(8000) # m.s-1
(22.544275306567194, 140.17120347383343)

>>> plane.compute_stall_speed(8000, C_L_max=1.5) # m.s-1
41.80281924283373

Here is another example with the ascension speed and slop at sea level for an almost empty tank:

>>> plane.m_fuel = 0  # kg
>>> plane.update_variables(True)
>>> plane.compute_max_ascension_speed(z=0) # m.s-1
32.89763560421959

>>> plane.compute_reference_speed(z=0) # m.s-1
34.523934888646956

>>> plane.compute_max_ascension_slope(z=0)
0.5695896796157822

Hypotheses
----------

All the computations performed in the model rely on some hypotheses. Here is the list of the main ones:

1. The flight is considered to be symmetric (no yaw angle). The aerodynamic forces therefore have only two components, the drag along the `x` axis and the lift along the `z` one.
2. The aerodynamic forces are functions of the surrounding flow velocity.
3. No wind effect is considered.
4. The thrust is parallel to the velocity vector and does not participate in the lift.

Please consult the documentation of each functions if needed for more details.