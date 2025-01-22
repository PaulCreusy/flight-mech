# Flight Mechanics Calculator

![PyPI - Version](https://img.shields.io/pypi/v/flight-mech)![PyPI - Downloads](https://img.shields.io/pypi/dm/flight-mech)![Pylint Badge](https://github.com/PaulCreusy/flight-mech/actions/workflows/pylint.yml/badge.svg)![Pytest Badge](https://github.com/PaulCreusy/flight-mech/actions/workflows/pytest.yml/badge.svg)

## License

This software has been developed by Paul Creusy and is shared under the MIT License.

## Getting started

### Installation

#### Pip installation

To install this module with pip, please use:

```bash
pip install flight-mech
```

#### Manual installation

For a manual installation, please clone the repository and install the required Python libraries using the command:

```bash
pip install -r requirements.txt
```

### Usage

This software includes a simple atmospheric model and a set of flight mechanics equations allowing to compute plane characteristics.

**Please note that all equations and variables are defined in the international unit system.**

The plane model allows to compute the following quantities:
- max glide ratio
- speed at specific angle of incidence and altitude
- drag
- lift
- thrust
- stall speed
- reference speed
- minimum descent gliding slope
- gliding speed
- maximum gliding time
- maximum gliding range
- authorized velocity interval at fixed thrust for flight at constant altitude
- thrust needed at fixed altitude and angle of incidence
- minimum thrust needed at fixed altitude
- speed at minimum thrust
- maximum flight altitude
- speed for maximum ascension speed
- ascension slope for a specific angle of incidence and altitude
- load factor in turn
- maximum range at fixed altitude
- maximum range at fixed speed
- endurance
- take off distance without friction
- take off distance with friction
- landing distance
- take off speed
- landing speed
- alpha and delta coefficient at a flight point

Additionally, the following graphs can be generated:
- polar graph
- thrust-speed graph
- power-speed graph

Some examples are provided in the `examples` folder (please note that they do not cover all the use cases) as well with a few plane models in the `plane_database` folder. 

Here is an overview of what the software can compute:

```python
# Load the plane
plane = Plane("cessna_172", "./plane_database")

# Compute the fmax and CL at fmax
print("C_L_f_max", plane.C_L_f_max)
print("fmax", plane.f_max)
```

```bash
>> C_L_f_max 0.7745966692414834
>> fmax 12.909944487358056
```

```python
# Compute the speed interval at 8000 meters
plane.m_fuel = 136.26  # kg
plane.update_variables(True)
print("reference speed at 8000m [m.s-1]", plane.compute_reference_speed(8000))
print("speed interval at 8000m [m.s-1]",
      plane.compute_velocity_interval_for_fixed_thrust(8000))
print("stall speed at 8000m [m.s-1]",
      plane.compute_stall_speed(8000, C_L_max=1.5))
```

```bash
>> reference speed at 8000m [m.s-1] 56.214394963985406
>> speed interval at 8000m [m.s-1] (22.544275306567194, 140.17120347383343)
>> stall speed at 8000m [m.s-1] 41.80281924283373
```

```python
# Compute the ascension speed and slope at sea level
plane.m_fuel = 0  # kg
plane.update_variables(True)
print("max ascension speed [m.s-1]", plane.compute_max_ascension_speed(z=0))
print("reference speed at 0m [m.s-1]", plane.compute_reference_speed(z=0))
print("max slope at 0m [%]", plane.compute_max_ascension_slope(z=0))
```

```bash
>> max ascension speed [m.s-1] 32.89763560421959
>> reference speed at 0m [m.s-1] 34.523934888646956
>> max slope at 0m [%] 0.5695896796157822
```

### Documentation

The documentation is available online [here](https://paulcreusy.github.io/flight-mech/).

Otherwise, if you decided to clone the repository, you can generate the documentation using the following commands:

```bash
cd docs
make html
```

And open the file `docs/_build/html/index.html` in your browser. 
