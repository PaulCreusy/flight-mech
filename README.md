# Flight Mechanics Calculator

## License

This software has been developed by Paul Creusy and is shared under the MIT License.

## Getting started

### Installation

To install this software, please clone the repository and install the required Python libraries using the command:

```bash
pip install -r requirements.txt
```

### Usage

This software includes a simple atmospheric model and a set of flight mechanics equations allowing to compute plane characteristics.

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