using CylindersBasedCameraResectioning.Geometry
using CylindersBasedCameraResectioning: sample_camera

@testset "issame_line" begin
	line₁ = Line([0, 0, 0], [1, 1, 1])
	line₂ = Line([1, 1, 1], [2, 2, 2])
	@test issame_line(line₁, line₂)
end

@testset "project_point_into_line" begin
	point = [1, 1, 1]
	line = Line([0, 0, 0], [0, 0, 1])
	projected_point = project_point_into_line(point, line)
	@test projected_point ≈ [0, 0, 1]
end

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

@testset "get_cylinder_contours" begin
	cylinder = Cylinder([-8.80334, 2.31298, 6.82885], 1, ([-7.30743, 2.1323, 10.412] - [-9.44826, 0.687196, 3.35754]))
	camera_center = sample_camera.position
	camera_matrix = sample_camera.matrix
	camera_matrix = camera_matrix ./ camera_matrix[3, 4]
	contour₁, contour₂ = get_cylinder_contours(cylinder, camera_center, camera_matrix)

	point₁₁ = [-10.3515, 0.769649, 3.61477]
	point₁₂ = [-8.18179, 2.18732, 10.663]
	point₂₁ = [-9.13133, 2.1864, 2.96598]
	point₂₂ = [-6.9397, 3.58753, 9.99881]

	projected_point₁₁ = camera_matrix * [point₁₁; 1]
	projected_point₁₂ = camera_matrix * [point₁₂; 1]
	projected_point₂₁ = camera_matrix * [point₂₁; 1]
	projected_point₂₂ = camera_matrix * [point₂₂; 1]

	projected_point₁₁ = projected_point₁₁ ./ projected_point₁₁[3]
	projected_point₁₂ = projected_point₁₂ ./ projected_point₁₂[3]
	projected_point₂₁ = projected_point₂₁ ./ projected_point₂₁[3]
	projected_point₂₂ = projected_point₂₂ ./ projected_point₂₂[3]

	projected_point₁₁ = projected_point₁₁[1:2]
	projected_point₁₂ = projected_point₁₂[1:2]
	projected_point₂₁ = projected_point₂₁[1:2]
	projected_point₂₂ = projected_point₂₂[1:2]

	target_contour₁ = Line(projected_point₁₁, projected_point₁₂ - projected_point₁₁)
	target_contour₂ = Line(projected_point₂₁, projected_point₂₂ - projected_point₂₁)

	@test issame_line(contour₁, target_contour₁, 1)
	@test issame_line(contour₂, target_contour₂, 1)
end