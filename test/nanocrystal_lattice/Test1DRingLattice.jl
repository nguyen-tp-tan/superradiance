using Test

project_path = "/home/lappy/Desktop/D_Baranov-L_Tan_Project/"
source_path = string(project_path, "Codes/superradiance/")

include(string(source_path, "NanocrystalLattice.jl"))

# function test_to_file(N, LNC, Lshell)
#   dir_name = string(@__DIR__, "/ring_1d")
#   if (~isdir(dir_name))
#     mkdir(dir_name)
#   end

#   file_name = string(dir_name,"N=",N,"LNC=",LNC,"nm_Lshell=",Lshell,"nm.txt")
#   open(file_name, "w") do file_address
#     writedlm(file_address, NanocrystalPositions1DRing(N, LNC, Lshell))
#   end
# end

if abspath(PROGRAM_FILE) == @__FILE__
  # test_to_file(10, 10.0, 1.5)
  # test_to_file(2, 10.0, 1.5)
  # test_to_file(1, 10.0, 1.5)
  # test_to_file(10, 10.0, 40.0)

  # Test the case where only 1 nanocrystal exists
  case1 = NanocrystalPositions1DRing(1, 10.0, 1.5)
  @test size(case1) == (3, 1)
  println("Size of first case = ", string(size(case1)))
  @test maximum(case1[3,:]) == 1.5 && minimum(case1[3,:]) == 1.5

  # Test the case where there are 2 nanocrystals along x-axis
  case2 = NanocrystalPositions1DRing(2, (4.0+2*1)*pi, 1.0)
  @test size(case2) == (3, 2)
  println("Size of second case = ", string(size(case2)))
  @test case2[1:2, 1] == [6.0, 0.0]
  @test case2[1:2, 2] == [-6.0, 0.0]
  @test maximum(case2[3,:]) == 1.0 && minimum(case2[3,:]) == 1.0

  # Test the case where there are 4 nanocrystals 
  # 2 along x-axis, 2 along y-axis
  case3 = NanocrystalPositions1DRing(4, (2.0+2*0.25)*pi, 0.0)
  @test size(case3) == (3, 4)
  println("Size of third case = ", string(size(case3)))
  @test case3[1:2, 1] == [5.0, 0.0]
  @test case3[1:2, 2] == [0.0, 5.0]
  @test case3[1:2, 3] == [-5.0, 0.0]
  @test case3[1:2, 4] == [0.0, -5.0]
  @test maximum(case3[3,:]) == 0.0 && minimum(case3[3,:]) == 0.0

  # Test the case with an even number of nanocrystals
  # The nanocrystal with index = N/2 + 1 is along x-axis
  N = 2*20
  case4 = NanocrystalPositions1DRing(N, 13.0, 0.75)
  @test size(case4) == (3, N)
  println("Size of fourth case = ", string(size(case4)))
  @test case4[2, 1] == 0.0
  @test case4[2, N ÷ 2 + 1] == 0.0
  @test maximum(case4[3, :]) == 0.75 && minimum(case4[3, :]) == 0.75

  # Test the case with multiple of 12, i.e. mult*12
  # The y-coordinate = 1/2 radius for the (mult+1) point
  # The x-coordinate = 1/2 radius for the (2mult+1) point
  mult = 10
  N = mult*12
  case5 = NanocrystalPositions1DRing(N, 4.0*pi, 0.5*pi)
  println(case5[2,mult+1])
  println(case5[1,2*mult+1])
  @test maximum(case5[2,mult+1] - 12*mult) < eps(Float64) # 1.0e-15
  @test maximum(case5[1,2*mult+1] - 12*mult) < eps(Float64) # 1.0e-15
  @test abs(maximum(case5[3,:])-0.5*pi) < eps(Float64) && abs(minimum(case5[3,:])-0.5*pi) < eps(Float64)
end
