using Test

project_path = "/home/lappy/Desktop/D_Baranov-L_Tan_Project/"
source_path = string(project_path, "Codes/superradiance/")

include(string(source_path, "NanocrystalLattice.jl"))
include(string(source_path, "ConstructHamiltonian.jl"))

function test_get_empty_sites()
  n_empty = 1
  n_link = 1
  nx = n_empty*2+n_link
  indices_of_suplatt =  NanocrystalIndexing(nx, nx, 1)
  empty_sites = get_empty_sites(n_empty, n_link, indices_of_suplatt)
  @test empty_sites == vec([3 7])

  n_empty = 1
  n_link = 2
  nx = n_empty*2+n_link
  indices_of_suplatt =  NanocrystalIndexing(nx, nx, 1)
  empty_sites = get_empty_sites(n_empty, n_link, indices_of_suplatt)
  @test empty_sites == vec([4 13])

  n_empty = 2
  n_link = 3
  nx = n_empty*2+n_link
  indices_of_suplatt =  NanocrystalIndexing(nx, nx, 1)
  empty_sites = get_empty_sites(n_empty, n_link, indices_of_suplatt)
  @test empty_sites == vec([6 7 13 14 36 37 43 44])
end

function count_zero_matrix_entries(n_empty, n_link, indices)
  n = (n_empty * 2 + n_link)^2
  empty_sites = get_empty_sites(n_empty, n_link, indices)

  number_of_zeros = 0
  for irow=1:n
    for icol=1:n
      is_row_empty = typeof(findfirst(i->i==irow, empty_sites)) != Nothing
      is_col_empty = typeof(findfirst(i->i==icol, empty_sites)) != Nothing
      if (is_row_empty || is_col_empty)
        number_of_zeros += 1
      end
    end
  end
  return number_of_zeros
end

function test_count_zeros_of_matrices()
  n_empty = 1
  n_link = 1
  nx = n_empty*2+n_link
  indices_of_suplatt =  NanocrystalIndexing(nx, nx, 1)
  numb_of_zeros = count_zero_matrix_entries(n_empty, n_link, indices_of_suplatt)
  @test numb_of_zeros == 32

  n_empty = 1
  n_link = 2
  nx = n_empty*2+n_link
  indices_of_suplatt =  NanocrystalIndexing(nx, nx, 1)
  numb_of_zeros = count_zero_matrix_entries(n_empty, n_link, indices_of_suplatt)
  @test numb_of_zeros == 60

  n_empty = 3
  n_link = 1
  nx = n_empty*2+n_link
  indices_of_suplatt =  NanocrystalIndexing(nx, nx, 1)
  numb_of_zeros = count_zero_matrix_entries(n_empty, n_link, indices_of_suplatt)
  @test numb_of_zeros == 49^2 - 31^2
end

function test_number_of_zeros_of_dumbbell_Hamiltonian()
  n_empty = 1
  n_link = 1
  nx = n_empty*2+n_link
  indices_of_suplatt =  NanocrystalIndexing(nx, nx, 1)
  kr = NanocrystalPositions(nx, nx, 1, 1.0, indices_of_suplatt)
  H_total = make_2d_dumbbell_Htotal(n_empty, n_link, 1.0, indices_of_suplatt, kr, 
                                    1.0, zeros(nx^2,1), radiative_coupling_full)
  zero_matrice_elements = findall(==(0), H_total)
  @test size(zero_matrice_elements)[1] = 32
  @test Set(zero_matrice_elements) == Set(
    [CartesianIndex(1, 3), CartesianIndex(1,7), CartesianIndex(3, 1), CartesianIndex(7,1),
     CartesianIndex(2, 3), CartesianIndex(2,7), CartesianIndex(3, 2), CartesianIndex(7,2),
     CartesianIndex(4, 3), CartesianIndex(4,7), CartesianIndex(3, 4), CartesianIndex(7,4),
     CartesianIndex(5, 3), CartesianIndex(5,7), CartesianIndex(3, 5), CartesianIndex(7,5),
     CartesianIndex(6, 3), CartesianIndex(6,7), CartesianIndex(3, 6), CartesianIndex(7,6),
     CartesianIndex(8, 3), CartesianIndex(8,7), CartesianIndex(3, 8), CartesianIndex(7,8),
     CartesianIndex(9, 3), CartesianIndex(9,7), CartesianIndex(3, 9), CartesianIndex(7,9),
     CartesianIndex(3, 3), CartesianIndex(3,7), CartesianIndex(7, 3), CartesianIndex(7,7)])

  n_empty = 1
  n_link = 2
  nx = n_empty*2+n_link
  indices_of_suplatt =  NanocrystalIndexing(nx, nx, 1)
  kr = NanocrystalPositions(nx, nx, 1, 1.0, indices_of_suplatt)
  H_total = make_2d_dumbbell_Htotal(n_empty, n_link, 1.0, indices_of_suplatt, kr, 
                                    1.0, zeros(nx^2,1), radiative_coupling_full)
  zero_matrice_elements = findall(==(0), H_total)
  @test size(zero_matrice_elements)[1] = 60

  n_empty = 3
  n_link = 1
  nx = n_empty*2+n_link
  indices_of_suplatt =  NanocrystalIndexing(nx, nx, 1)
  kr = NanocrystalPositions(nx, nx, 1, 1.0, indices_of_suplatt)
  H_total = make_2d_dumbbell_Htotal(n_empty, n_link, 1.0, indices_of_suplatt, kr, 
                                    1.0, zeros(nx^2,1), radiative_coupling_full)
  zero_matrice_elements = findall(==(0), H_total)
  @test size(zero_matrice_elements)[1] = 49^2 - 31^2
end

if abspath(PROGRAM_FILE) == @__FILE__
  test_get_empty_sites()
  test_count_zeros_of_matrices()
  test_number_of_zeros_of_dumbbell_Hamiltonian
end