# This file contains the helper/utility functions

# ------------------------------------------

"""
Screening factor for a dielectric sphere inside a dielectric medium
```julia
julia> spherical_screening_factor(1.2, 5.1)
```
## INPUT
`eps_out`: dielectric constant of the external medium
`eps_in`: dielectric constant of the inner material
"""
function spherical_screening_factor(eps_out, eps_in)
  return 3.0 / (2.0 + eps_in/eps_out)
end


"""
Ratio of the screening factors when the dielectric constant of the external medium is changed
```julia
julia> spherical_screening_factor(1.2, 2.0, 5.1)
```
## INPUT
`eps_o1`: dielectric constant of the first external medium
`eps_o2`: dielectric constant of the second external medium
`eps_in`: dielectric constant of the inner material
"""
function rate_ratio_spherical_by_varying_external_dielectric(eps_o1, eps_o2, eps_in)
  f_eps_1 = spherical_screening_factor(eps_o1, eps_in)
  f_eps_2 = spherical_screening_factor(eps_o2, eps_in)
  return sqrt(eps_o1/eps_o2) * (f_eps_1/f_eps_2)^2
end


"""
Ratio of the screening factors when the dielectric constant of the external medium is changed
```julia
julia> spherical_screening_factor(1.2, 4.5, 5.0)
```
## INPUT
`eps_out`: dielectric constant of the external medium
`eps_i1`: dielectric constant of the first inner material
`eps_i2`: dielectric constant of the second inner material
"""
function rate_ratio_spherical_by_varying_inner_dielectric(eps_out, eps_i1, eps_i2)
  f_eps_1 = spherical_screening_factor(eps_out, eps_i1)
  f_eps_2 = spherical_screening_factor(eps_out, eps_i2)
  return (f_eps_1/f_eps_2)^2
end

# ------------------------------------------------------------------------------------
