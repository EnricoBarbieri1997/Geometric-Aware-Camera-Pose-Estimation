const TESTS = true

using Test

@testset "Cylinder camera resectioning" begin
    @testset "Utils" begin
        include("utils_tests.jl")
    end
    @testset "Geometry" begin
        include("geometry_tests.jl")
    end
    @testset "Camera" begin
        include("camera_tests.jl")
    end
    @testset "Data representation" begin
        include("data_representation_tests.jl")
    end
    @testset "IO" begin
        include("io_tests.jl")
    end
    @testset "Equation Systems" begin
        include("equation_system_tests.jl")
    end
end
