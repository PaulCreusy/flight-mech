# Flight-mech

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

### Documentation

The documentation is available online [here](https://flight-mech.creusy.fr).

Otherwise, if you decided to clone the repository, you can generate the documentation using the following commands:

```bash
cd docs
make html
```

And open the file `docs/_build/html/index.html` in your browser. 

### Functionalities

This software includes various modules to build a numerical plane model and compute its characteristics. The modules implemented are the following:

- `aerodynamics` : contains functions to compute quantities in the boundary layer of a fluid flow.
- `airfoil` : allows to define the geometry of an airfoil and compute the lift and moment coefficients.
- `atmosphere` : defines several atmosphere models to compute density, temperature, pressure and other quantities.
- `environment` : defines environmental quantities such as gravity.
- `fuel` : defines several types of broadly used fuels in aeronautics.
- `gas` : defines several gas models to compute their properties.
- `maneuver` : allows to define maneuvers for a plane model.
- `motor` : defines an electric motor model.
- `plane` : allows to define a numerical plane model, binding the previous modules, to compute its flight characteristics.
- `turbine` : allows to define several types of turbine to compute their thrust and consumption at various operating conditions.
- `wing` : allows to define the geometry of a wing and compute the lift and drag coefficients.
