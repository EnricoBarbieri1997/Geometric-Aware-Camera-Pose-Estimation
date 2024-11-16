module EquationSystems
	export stack_intrinsic_rotation_conic_parameters, build_intrinsic_rotation_conic_system, build_intrinsic_rotation_translation_conic_system

	using ..Camera: build_camera_matrix
	using ..Space: build_rotation_matrix

	using HomotopyContinuation

	function stack_intrinsic_rotation_conic_parameters(lines, points_at_infinity)
		return vcat(vec(lines), vec(points_at_infinity))
	end

	function build_intrinsic_rotation_conic_system(lines_values::Matrix{<:Number}, points_at_infinity_values::Matrix{<:Number})
		lines_count = size(lines_values)[1]
		@var x y z f
		@var lines[1:lines_count, 1:3] points_at_infinity[1:lines_count, 1:3]
		Rₚ = build_rotation_matrix(x, y, z)

		system_to_solve = []
		for line_index in 1:lines_count
			equation = lines[line_index, :]' * [
					f 0 0;
					0 f 0;
					0 0 1
			] * Rₚ * points_at_infinity[line_index, :]
			push!(system_to_solve, equation)
		end
		parameters = stack_intrinsic_rotation_conic_parameters(lines, points_at_infinity)
		return System(system_to_solve, variables = [x, y, z, f], parameters = parameters)
	end

	function build_intrinsic_rotation_translation_conic_system(intrinsic, rotation, lines_values, dual_quadic_values)
		lines_count = size(lines_values)[1]
		@var tx ty tz
		@var lines[1:lines_count, 1:3] dual_quadrics[1:lines_count, 1:4, 1:4]
		P = build_camera_matrix(intrinsic, rotation, [tx; ty; tz])

		system_to_solve = []
		for i in 1:lines_count
				equation = lines[i, :]' * P * dual_quadrics[i, :, :] * P' * lines[i, :]
				push!(system_to_solve, equation)
		end
		parameters = stack_intrinsic_rotation_conic_parameters(lines, dual_quadrics)
		return System(system_to_solve, variables=[tx, ty, tz], parameters=parameters)
	end
end