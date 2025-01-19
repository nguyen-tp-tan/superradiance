using Test

include("../../Helper.jl")

function test_spherical_screening_factor()
  eps_in = 2.0
  eps_out = eps_in
  @test spherical_screening_factor(eps_out, eps_in) == 1.0

  eps_in = 4.0
  eps_out = 1.0
  @test spherical_screening_factor(eps_out, eps_in) == 0.5
end

function test_rate_ratio_by_varying_external_dielectric()
  eps_in = 4.0
  eps_out1 = 2.0
  eps_out2 = eps_out1
  @test rate_ratio_spherical_by_varying_external_dielectric(eps_out1, eps_out2, eps_in) == 1.0

  eps_in = 4.0
  eps_out1 = 4.0
  eps_out2 = 1.0
  @test rate_ratio_spherical_by_varying_external_dielectric(eps_out1, eps_out2, eps_in) == 8.0

  eps_in = 6.0
  eps_out1 = 3.0
  eps_out2 = 1.5
  @test rate_ratio_spherical_by_varying_external_dielectric(eps_out1, eps_out2, eps_in) == sqrt(2)*1.5^2
end

function test_rate_ratio_by_varying_inner_dielectric()
  eps_in1 = 5.0
  eps_in2 = eps_in1
  eps_out = 1.5
  @test rate_ratio_spherical_by_varying_inner_dielectric(eps_out, eps_in1, eps_in2) == 1.0

  eps_in1 = 6.0
  eps_in2 = 3.0
  eps_out = 1.5
  @test rate_ratio_spherical_by_varying_inner_dielectric(eps_out, eps_in1, eps_in2) == (2.0/3.0)^2
end

if abspath(PROGRAM_FILE) == @__FILE__
  test_spherical_screening_factor()
  test_rate_ratio_by_varying_external_dielectric()
  test_rate_ratio_by_varying_inner_dielectric()
end
