using CylindersBasedCameraResectioning
using Test

@testset "Cylinder camera resectioning" begin
	@testset "Utils" begin
		include("utils_tests.jl")
	end
	@testset "Geometry" begin
		include("geometry_tests.jl")
	end
end