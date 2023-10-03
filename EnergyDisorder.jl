# This file contains functions required for studying energy disorder

using Random
using Distributions
    
# ------------------------------------------
"""
Function that calculates the standard deviation of energy distribution
taking kinetic energy, Coulomb energy and size dispersion as input
```julia
julia> EneryWidth(0.5, -0.3, 0.1, 1)
```
## INPUT
`Ekin`: kinetic energy (energy unit)
`ECoul`: Coulomb energy (energy unit)
`size_stand_dev`: standard deviation of the size dispersion (%)
`max_deriv`: highest order of derivative to compute the energy spread
"""
function EnergyWidth(Ekin, ECoul, size_stand_dev, max_deriv)
    result = 0.0
    for i_deriv=1:max_deriv
        result = result + (-1*size_stand_dev)^i_deriv * ((i_deriv+1)*Ekin + ECoul)
    end
    return abs(result)
end
    
# ------------------------------------------

"""
Generate an array following gaussian distribution
``` julia
julia> gen_gaussian_rand_array(0.5, 100)
```
`mean_value`: mean value of the Gaussian distribution
`standard_deviation`: standard deviation for the Gaussian distribution
`numb_element`: number of elements (size) of output array
"""
function gen_gaussian_rand_array(mean_value, standard_deviation, numb_element)
    Gauss_obj = Normal(mean_value, standard_deviation)
    rand_array = rand(Gauss_obj, numb_element)
    return rand_array
end
