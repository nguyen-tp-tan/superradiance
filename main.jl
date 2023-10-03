# initialization of packages and functions

using LinearAlgebra
using DelimitedFiles

project_path = "/home/lappy/Desktop/D_Baranov-L_Tan_Project"
source_path = joinpath(project_path, "Codes/superradiance")
vec_soure_files = ["Constants.jl", "ConstructHamiltonian.jl", 
"NanocrystalLattice.jl", "EnergyDisorder.jl"]
for i_file = 1:size(vec_soure_files)[1]
  include(joinpath(source_path, vec_soure_files[i_file]))
end

# ---------------------------------------------------------------------------------

# Inputs for radiative lifetime, energy, dielectric values, size
# tau_r: radiative lifetime, in nanosecond
# omega0: bright single-exciton energy, in eV
# eps_opt: optical dielectric constant
# LNC: average edge length of a nanocrystal (nm)
# Lshell: average length of ligand shell (nm)
tau_r = 1.0 
omega0 = 2.636
eps_opt = 5.0
LNC = 3.8
Lshell = 0.35

# gamma_r: radiative decay rate, in eV
# E_diag: diagonal energy values including imaginary part for decay, in eV
# kLNC2NC: center-to-center distance between adjacent nanocrystals
# prefac: prefactor for the radiative coupling, = 1/2 * gamma_r
gamma_r, E_diag, kLNC2NC = initialize_parameters(tau_r, omega0, eps_opt, LNC, Lshell)
prefac = 0.5*gamma_r

dirname_LNC = joinpath(pwd(), string("LNC=", LNC, "nm_Lshell=", Lshell, "nm"))
if ~isdir(dirname_LNC); mkdir(dirname_LNC) end 

# Inputs for energy disorder
# delta_E: standard deviation (eV) for the Gaussian distribution of energies of NCs
delta_E = 0.0242 

# nz_min: minimum value of nz
# nz_max: maximum value of nz
# nz_step: step size for iterating from nz_min to nz_max
nz_min = 1 
nz_max = 40 
nz_step = 1 

# --------------------------------------------------------------------------------

for nz=nz_min:nz_step:nz_max 
  nx = 1; ny = nz;
  nNC = nx*ny*nz

  # index_of_NC: indices (ix,iy,iz) for all nanocrystals
  # kr_of_NC: k*position (r_x,r_y,r_z) for all nanocrystals
  index_of_NC =  NanocrystalIndexing(nx, ny, nz)
  kr_of_NC = NanocrystalPositions(nx, ny, nz, kLNC2NC, index_of_NC)

  # vec_delta_E: vector of delta_E values for each diagonal block
  vec_delta_E = gen_gaussian_rand_array(0, delta_E, nNC)
  
  # Construct the Hamiltonian with radiative coupling
  H_total = make_radiative_Htotal(nNC, prefac, kr_of_NC, E_diag, vec_delta_E)
  
  # ---------------------------------------------------------------------------------
  
  # solve the total Hamiltonian, store results for each superlattice in .txt file
  # eig_H_total: eigen values of the total Hamiltonian H_total
  # Energy_k: energies, real part of eigen states of total Hamiltonian H_total
  # Gamma_k: all decay rates (positive values), (-1) * imaginary part of eigen values 
  eig_H_total = eigvals(H_total)
  Energy_k = real(eig_H_total)
  Gamma_k = -imag(eig_H_total)

  fname_nNC = joinpath(dirname_LNC, string("SupLatt_", nx, "x", ny, "x", nz, ".txt"))
  fout_id = open(fname_nNC,"w")
  println(fout_id,"## E_k \t Gamma_k:")
  for ielement=1:3*nNC
    println(fout_id, Energy_k[ielement], "\t", Gamma_k[ielement])
  end
  close(fout_id)

  H_total = Nothing
end
