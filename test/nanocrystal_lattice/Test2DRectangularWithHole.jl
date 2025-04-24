using Test

project_path = "/home/lappy/Desktop/D_Baranov-L_Tan_Project/"
source_path = joinpath(project_path, "Codes/superradiance/")

include(joinpath(source_path, "NanocrystalLattice.jl"))

if abspath(PROGRAM_FILE) == @__FILE__
  # empty matrix
  mat_nul = NanocrystalPositions2DRectangularWithHole(20, 20, 0.5, 0.25, 0.0, 1.01, 0.90)
  @test mat_nul == zeros(3, 0)

  # full rectangle
  mat_full = NanocrystalPositions2DRectangularWithHole(4, 5, 0.5, 0.25, 0.0, 0.0, 0.0)
  @test size(mat_full) == (3, 20)
  @test mat_full[1:3,3] == [1.5; 0.25; 0.0]
  @test mat_full[1:3,14] == [1.0; 1.00; 0.0]

  # rectangular strips along x
  mat_xstrips = NanocrystalPositions2DRectangularWithHole(10, 3, 0.25, 0.5, 0.0, 1.0, 0.3)
  @test size(mat_xstrips) == (3, 20)
  @test mat_xstrips[1:3,1] == [0.25; 0.5; 0.0]
  @test mat_xstrips[1:3,7] == [7*0.25; 0.5; 0.0]
  @test mat_xstrips[1:3,12] == [2*0.25; 1.5; 0.0]

  # rectangular strips along y
  mat_ystrips = NanocrystalPositions2DRectangularWithHole(10, 3, 0.25, 0.5, 1.2, 0.3, 1.0)
  @test size(mat_ystrips) == (3, 24)
  @test mat_ystrips[1:3,8] == [2.5; 0.5; 1.2]
  @test mat_ystrips[1:3,11] == [0.75; 1.0; 1.2]
  @test mat_ystrips[1:3,20] == [1.0; 1.5; 1.2]

  # rectangular lattice with a hole
  mat_hole = NanocrystalPositions2DRectangularWithHole(5, 4, 1.00, 2.0, 1.0, 0.3, 0.4)
  @test size(mat_hole) == (3, 18)
  @test mat_hole[1:3, 1] == [1.0; 2.0; 1.0]
  @test mat_hole[1:3, 5] == [5.0; 2.0; 1.0]
  @test mat_hole[1:3, 7:8] == [2.0 4.0; 4.0 4.0; 1.0 1.0]
  @test mat_hole[1:3, 11:12] == [2.0 4.0; 6.0 6.0; 1.0 1.0]
  @test mat_hole[1:3, 18] == [5.0; 8.0; 1.0]
end
