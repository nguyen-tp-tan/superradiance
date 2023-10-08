# This file contains all the functions needed for the radiative coupling term in the total Hamiltonian.
# The details of each function can be found in the comment of each function below.

using GSL

# ------------------------------------------

"""
Function that computes the following matrix
```
[(ehat_alpha.ehat_beta) - 3(ehat_alpha.rhat)(ehat_beta.rhat)]
```
```julia
julia> DipSecondTerm([1.0 0.0 1.0])
```
## Input:
`rhat` : normalized 3-tuple = unit vector of r
"""
function DipSecondTerm(rhat)
  size_of_r = size(rhat)
  if (minimum(size_of_r)!=1 || maximum(size_of_r)!=3)
    result = 0.0
  else
    if (size_of_r[2]==1) rhat = reshape(rhat,1,3) end
    result = Matrix{Float64}(I,3,3) - 3.0 * (rhat' * rhat)
  end
  return result 
end

# ------------------------------------------

"""
Function that computes the full radiative coupling with J_{i,j}

### INPUTS:
`norm_kr`: norm of k*r, k=wave number of the transition, r=position vector

`rhat`: 3x1 unit vector along the direction of r
"""
function radiative_coupling_full(norm_kr, rhat)
  y0_j0 = sf_bessel_y0(norm_kr) - im * sf_bessel_j0(norm_kr) 
  y2_j2 = 0.5 * ( im * sf_bessel_j2(norm_kr) - sf_bessel_y2(norm_kr) )
  result = y0_j0 .* Matrix{Float64}(I,3,3) .+ y2_j2 .* DipSecondTerm(rhat)
  return result
end

# ------------------------------------------

"""
Function that computes the off-diagonal dipole-only y2 terms of Omega_{i,j}

### INPUTS:

`norm_kr`: norm of k*r, k=wave number of the transition, r=position vector

`rhat`: 3x1 unit vector along the direction of r
"""
function dipole_coupling_y2(norm_kr, rhat)
  y2 = - 0.5 * sf_bessel_y2(norm_kr)
  result = y2 .* DipSecondTerm(rhat)
  return result
end

# ------------------------------------------

"""
Function that computes the off-diagonal dipole-mediated j2 and y2 terms of Omega_{i,j}

### INPUTS:

`norm_kr`: norm of k*r, k=wave number of the transition, r=position vector

`rhat`: 3x1 unit vector along the direction of r
"""
function dipole_coupling_j2y2(norm_kr, rhat)
  y2_j2 = 0.5 * ( im * sf_bessel_j2(norm_kr) - sf_bessel_y2(norm_kr) )
  result = y2_j2 .* DipSecondTerm(rhat) 
  return result
end

# ------------------------------------------

"""
Function that computes the off-diagonal dipole-mediated coupling 
```julia
julia> RadiativeCoupling(radiative_coupling_full, 0.5, [1/2, 1/2, 0])
```
### INPUTS:
`coupling_func`: function name for which coupling function to be used

`prefac`: prefactor = hbar*gamma_r/2; gamma_r = radiative decay rate

`kr`: k*r, k=wave number of the transition, r=position vector
"""
function RadiativeCoupling(coupling_func::Function, prefac, kr)
  size_kr = size(kr)
  if (minimum(size_kr)!=1 || maximum(size_kr)!=3)
    result = 0.0
  else
    if (size_kr[2]==1) kr = reshape(kr,1,3) end
    norm_kr = sqrt(dot(kr,kr))
    if (norm_kr > eps(Float64))
      rhat = kr ./ norm_kr
      result = prefac .* coupling_func(norm_kr, rhat)
    else
      result = Inf
    end
  end
  return result
end

