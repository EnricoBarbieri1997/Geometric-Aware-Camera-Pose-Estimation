using CylindersBasedCameraResectioning.Utils

@testset "almostequal" begin
	@test almostequal(1e-7, 0) == true
	@test almostequal(1e-5, 0) == false
end