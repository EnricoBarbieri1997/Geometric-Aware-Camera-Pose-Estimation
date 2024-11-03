module EquationSystems
	using ..Camera: build_camera_matrix
	using ..Space: build_rotation_matrix

	function build_intrinsic_rotation_conic_system(lines_values, points_at_infinity_values)
		lines_count = length(lines_values)
		@var x y z f
		@var lines[1:lines_count, 1:3] points_at_infinity[1:lines_count, 1:3]
		Rₚ = build_rotation_matrix(x, y, z)

		system_to_solve = []
		for (line_index, _) in lines
				equation = lines[line_index]' * [
						f 0 0;
						0 f 0;
						0 0 1
				] * Rₚ * points_at_infinity[line_index]
				push!(system_to_solve, equation)
		end
		return System(system_to_solve, variables = [x, y, z, f], parameters = [lines, points_at_infinity])
	end

	function build_intrinsic_rotation_translation_conic_system(intrinsic, rotation, lines_values, dual_quadic_values)
		lines_count = length(lines_values)
		@var tx ty tz
		@var lines[1:lines_count, 1:3] dual_quadrics[1:lines_count, 1:4, 1:4]
		P = build_camera_matrix(intrinsic, rotation, [tx; ty; tz])

		system_to_solve = []
		for i in 1:lines_count
				equation = lines[i]' * P * dual_quadrics[i] * P' * lines[i]
				push!(system_to_solve, equation)
		end

		return System(system_to_solve, variables=[tx, ty, tz], parameters=[lines, dual_quadrics])
	end
end