using Test

project_path = "/home/lappy/Desktop/D_Baranov-L_Tan_Project/"
source_path = string(project_path, "Codes/superradiance/")

include(string(source_path, "NanocrystalLattice.jl"))


if abspath(PROGRAM_FILE) == @__FILE__
  # Subsequent cases has unit distance a = 1 between adjacent points
  # Test the case where n2 = 1, i.e. only points along x-axis
  n1 = 10
  array1 = NanocrystalPositions2DHexagonal(1.0, n1, 1, 0.0)
  @test size(array1) == (3, 10)
  @test array1[1,:] == [1.0*i for i in 0:(n1-1)]
  @test maximum(array1[2,:]) == 0.0 && minimum(array1[2, :]) == 0.0
  @test maximum(array1[3,:]) == 0.0 && minimum(array1[3, :]) == 0.0 

  # Test the case where n1 = 1, i.e. only zig-zag points along y-axis
  n2 = 12
  array2 = NanocrystalPositions2DHexagonal(1.0, 1, n2, 1.2)
  @test size(array2) == (3, 12)
  @test array2[1,:] == [0.5*(i%2) for i=0:(n2-1)]
  @test array2[2,:] == sqrt(3) * [i//2 for i=0:(n2-1)]
  @test maximum(array2[3,:]) == 1.2 && minimum(array2[3, :]) == 1.2

  # Test the case where the hexagonal lattice confined 
  # within a rectangular boundary
  n1 = 8
  n2 = 11
  zValue = 2.05
  array3 = NanocrystalPositions2DHexagonal(1.0, n1, n2, zValue)
  @test size(array3) == (3, 88)
  @test maximum(array3[3,:]) == zValue && minimum(array3[3,:]) == zValue
  for irow=1:n2
    yValue = sqrt(3.0) * ((irow-1)//2)
    startIdx = (irow-1)*n1 + 1
    endIdx = (irow-1)*n1 + n1
    @test maximum(array3[2,startIdx:endIdx]) == yValue && minimum(array3[2,startIdx:endIdx]) == yValue
    array3Xslice = 0.5*mod((irow-1), 2) .+ [icol-1 for icol=1:n1]
    @test array3[1,startIdx:endIdx] == circshift(array3Xslice, floor(irow/2))
  end
end