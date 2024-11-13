using CylindersBasedCameraResectioning.Geometry

@testset "project_point_into_plane" begin
	point = [1, 1, 1]
	plane = Plane([0, 0, 0], [0, 0, 1])
	projected_point = project_point_into_plane(point, plane)
	@test projected_point ≈ [1, 1, 0]
end

@testset "get_tangentpoints_circle_point" begin
	circle = Circle([1, 3, 2], 4, [0, -2, -1])
	point = [-6, 3, 2]
	correctpoint₁ = [-1.29, 4.47, -0.94]
	correctpoint₂ = [-1.29, 1.53, 4.94]
	tangentpoint₁, tangentpoint₂ = get_tangentpoints_circle_point(circle, point)
	tangentpoint₁ = round.(tangentpoint₁; digits=2)
	tangentpoint₂ = round.(tangentpoint₂; digits=2)

	@test (tangentpoint₁ == correctpoint₁ && tangentpoint₂ == correctpoint₂) ||
		(tangentpoint₂ == correctpoint₁ && tangentpoint₁ == correctpoint₂)
end