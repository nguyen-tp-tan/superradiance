# initialization of packages and functions

using LinearAlgebra
using DelimitedFiles

include("ConstructHamiltonian.jl")
include("EnergyDisorder.jl")
include("NanocrystalLattice.jl")
include("RadiativeDipoleCoupling.jl")

# ------------------------------------------

"""
Function main to be executed only when the current script is executed

Step 1: provide all basic input parameters, initialize the remaining parameters
        
        tau_r: radiative lifetime (nanoseconds)

        omega0: bright single-exciton energy (eV)

        eps_opt: optical dielectric constant

        LNC: average edge length of a nanocrystal (nm)

        Lshell: average length of ligand shell (nm)

        gamma_r: radiative decay rate (eV)

        E_diag: diagonal values including imaginary part representing radiative decay (eV)

        kLNC2NC: center-to-center separation between adjacent nanocrystals

        prefac: prefactor of the radiative coupling, = 1/2 * gamma_r

        delta_E: standard deviation (eV) of the Gaussian distribution of energies of NCs
        
        nz_min: minimum value of nz
        
        nz_max: maximum value of nz
        
        nz_step: step size of iteration from nz_min to nz_max

Step 2: loop over superlattice size where there are 2 major parts

        2.1 Compute total Hamiltonian within radiative coupling model
        
            index_of_NC: indices (ix,iy,iz) of all nanocrystals

            kr_of_NC: k*position (r_x,r_y,r_z) of all nanocrystals

            vec_delta_E: vector of delta_E values for each diagonal block

            H_total: radiative Hamiltonian of nanocrystal superlattice system

        2.2 Compute eigen values/eigen vectors and write output to file for each case

            eig_H_total: eigen values of the total Hamiltonian H_total

            Energy_k: energies, real part of eigen states of total Hamiltonian H_total

            Gamma_k: all decay rates (positive values), (-1) * imaginary part of eigen values 

"""
function main()
  tau_r = 1.0 
  omega0 = 2.636
  eps_opt = 5.0
  LNC = 3.8
  Lshell = 0.35

  gamma_r, E_diag, kLNC2NC = initialize_parameters(tau_r, omega0, eps_opt, LNC, Lshell)
  prefac = 0.5*gamma_r

  dirname_LNC = joinpath(pwd(), string("LNC=", LNC, "nm_Lshell=", Lshell, "nm"))
  if ~isdir(dirname_LNC); mkdir(dirname_LNC) end 

  delta_E = 0.0242 

  nz_min = 1 
  nz_max = 25 
  nz_step = 1 

  # ------------------------------------------

  for nz=nz_min:nz_step:nz_max 
    nx = 1; ny = nz;
    nNC = nx*ny*nz

    index_of_NC =  NanocrystalIndexing(nx, ny, nz)
    kr_of_NC = NanocrystalPositions(nx, ny, nz, kLNC2NC, index_of_NC)

    vec_delta_E = gen_gaussian_rand_array(0, delta_E, nNC)

    H_total = make_radiative_Htotal(nNC, prefac, kr_of_NC, E_diag, vec_delta_E, radiative_coupling_full)
    
    # ------------------------------------------

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
end

# -------------------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
