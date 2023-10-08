# This file contains functions used in the construction of the total Hamiltonian

include("Constants.jl")
include("RadiativeDipoleCoupling.jl")

# ------------------------------------------

"""
Function to initialize important parameters: 
single-nanocrystal decay rate 'gamma_r', complex diagonal "energy" 'E_diag', 
center-to-center distance 'kd' between adjacent nanocrystals
from inputs:
lifetime, real resonant energy, optical dielectric constant, nanocrystal size, shell length
```julia
julia> initialize_parameters(1.0, 2.636, 5.0 3.8, 0.35)
```
## INPUT
`tau_r`: radiative lifetime (in nanosecond)

`omega`: resonant energy for nanocrystal of size 'LNC'

`eps_opt`: optical dielectric constant

`LNC`: nanocrystal size

`Lshell`: shell length (of ligand, encapsulation, etc.)
"""
function initialize_parameters(tau_r, omega, eps_opt, LNC, Lshell)
    gamma_r = (1/tau_r)*fnstoeV
    E_diag = omega - 0.5*im*(gamma_r)
    LNC2NC = LNC + 2*Lshell
    k_resonant = (omega/EautoeV)*sqrt(eps_opt)/(c_0au * Lautonm)
    kd = k_resonant*LNC2NC
    return gamma_r, E_diag, kd
end

# ------------------------------------------

"""
Function that constructs the total Hamiltonian containing 
the energy diagonal blocks and the off-diagonal radiative coupling term.
```julia
julia> make_radiative_Htotal(1000, 0.5, 2.3, zeros(1000, 1))
```
## INPUT
`n`: total number of nanocrystals

`prefac`: multiplicative pre-factor

`kr`: 3xn array of k_0 * r of the nanocrystals

`E_resonant`: the resonant energy for the diagonal blocks

`vec_delE`: vector of random energy variation, 
            = 0 if no energy disorder, != 0 otherwise
"""
function make_radiative_Htotal(n, prefac, kr, E_resonant, vec_delE, coupling_func::Function)
    Mat_I3x3 = Matrix{Complex{Float64}}(I, 3, 3)
    H_diag_blk = Mat_I3x3*E_resonant
    H_total = zeros(Complex{Float64}, 3*n, 3*n)
    for irow=1:n
      for icol=irow:n
        start_row = 3*(irow-1)+1; end_row = 3*irow
        start_col = 3*(icol-1)+1; end_col = 3*icol
        if (irow==icol)
          H_total[start_row:end_row,start_col:end_col] = H_diag_blk + Mat_I3x3*vec_delE[irow]
        else
          current_kr = kr[1:3,irow] - kr[1:3,icol]
          current_kr = reshape(current_kr, 3, 1)
          H_total[start_row:end_row,start_col:end_col] = RadiativeCoupling(coupling_func, prefac, current_kr)
          H_total[start_col:end_col,start_row:end_row] = H_total[start_row:end_row,start_col:end_col]
        end
      end
    end
    return H_total
end
