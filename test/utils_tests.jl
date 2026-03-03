using CylindersBasedCameraResectioning.Utils

using JSON

function parse_matrix(matrix_json)
	return Float64.(hcat(matrix_json...)')
end

function normalize_line(line)
	return line ./ line[3]
end

@testset "almostequal" begin
	@test almostequal(1e-7, 0) == true
	@test almostequal(1e-5, 0) == false
end

@testset "get_view_CC_VP_P" begin
	scene_path = joinpath(@__DIR__, "..", "assets", "test_scenes", "roller_coaster", "scene.json")
	views_path = joinpath(@__DIR__, "..", "assets", "test_scenes", "roller_coaster", "views.json")

	scene = open(scene_path) do io
		JSON.parse(io)
	end
	views = open(views_path) do io
		JSON.parse(io)
	end

	P = parse_matrix(scene["cameras"][1]["matrix"])
	view_lines = views[1]["lines"]

	for (index, cylinder) in enumerate(scene["cylinders"])
		CC = parse_matrix(cylinder["dual_matrix"])
		VP = Float64.(cylinder["singular_point"])
		ll1, ll2 = get_view_CC_VP_P(CC, VP, P)

		line1 = normalize_line(vec(ll1[:, 1]))
		line2 = normalize_line(vec(ll2[:, 1]))

		expected_lines = map(x -> normalize_line(Float64.(x)), view_lines["cylinder_$(index - 1)"])
		expected1 = expected_lines[1]
		expected2 = expected_lines[2]

		matches_direct = isapprox(line1, expected1; atol=1e-5) && isapprox(line2, expected2; atol=1e-5)
		matches_swapped = isapprox(line1, expected2; atol=1e-5) && isapprox(line2, expected1; atol=1e-5)
		@test matches_direct || matches_swapped
	end
end
