# Superradiance

Superradiance is a physics software package written in Julia programming language for studying superradiance in systems of nanocrystal superlattice.

## Installation

The package can be download directly from the Github Repo.
```bash
git clone git@github.com:nguyen-tp-tan/superradiance.git
```

## Usage

The package implements the radiative Hamiltonian as defined by, for instance, [Grad et al.] (http://dx.doi.org/10.1103/PhysRevA.37.3835).

Important functionalities are contained in:

    + "ConstructHamiltonian.jl": functions used to construct the total Hamiltonian of multi-nanocrystal superlattice

    + "EnergyDisorder.jl": functions used to describe parameters for the disorder system

    + "NanocrystalLattice.jl": functions for describing superlattice of nanocrystals

    + "RadiativeDipoleCoupling.jl": functions implementing dipole-mediated radiative coupling

Detailed documentation can be found in the respective files.

Run an example script, which is applied to a specific 2D system, by executing in Linux terminal 
```bash
julia main.jl
```

## Contribution

Pull requests are more than welcome. It would be nice to also open an issue to discuss what you would like to change.

## License
[GPLv3] (https://choosealicense.com/licenses/gpl-3.0/)
