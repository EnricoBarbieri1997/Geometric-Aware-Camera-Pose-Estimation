module EquationSystems
	export stack_homotopy_parameters, build_intrinsic_rotation_conic_system, build_intrinsic_rotation_translation_conic_system, build_camera_matrix_conic_system

	using ..Camera: build_intrinsic_matrix, build_camera_matrix
	using ..Space: build_rotation_matrix

	using HomotopyContinuation

	function stack_homotopy_parameters(parameters...)
		stacked_parameters = []
		for parameter in parameters
			stacked_parameters = vcat(stacked_parameters, vec(parameter))
		end
		return stacked_parameters
	end

	function build_intrinsic_rotation_conic_system(lines_values::Matrix{<:Number})
		lines_count = size(lines_values)[1]
		@var x y z f scale
		@var lines[1:lines_count, 1:3] points_at_infinity[1:lines_count, 1:3]
		Rₚ = build_rotation_matrix(x, y, z)
		intrinsic_topleft = build_intrinsic_matrix(f)[1:3, 1:3]

		system_to_solve = []
		for line_index in 1:lines_count
			equation = lines[line_index, :]' * scale * intrinsic_topleft * Rₚ * points_at_infinity[line_index, :]
			push!(system_to_solve, equation)
		end
		push!(system_to_solve, scale - 1)
		parameters::Vector{HomotopyContinuation.ModelKit.Variable} = stack_homotopy_parameters(lines, points_at_infinity)
		return System(system_to_solve, variables = [x, y, z, f, scale], parameters = parameters)
	end

	function build_intrinsic_rotation_translation_conic_system(intrinsic, rotation, lines_values)
		lines_count = size(lines_values)[1]
		@var tx ty tz scale
		@var lines[1:lines_count, 1:3] dual_quadrics[1:lines_count, 1:4, 1:4]
		P = build_camera_matrix(intrinsic, rotation, [tx; ty; tz])

		system_to_solve = []
		for i in 1:lines_count
				equation = lines[i, :]' * scale * P * dual_quadrics[i, :, :] * P' * lines[i, :]
				push!(system_to_solve, equation)
		end
		push!(system_to_solve, scale - 1)
		parameters::Vector{HomotopyContinuation.ModelKit.Variable} = stack_homotopy_parameters(lines, dual_quadrics)
		return System(system_to_solve, variables=[tx, ty, tz, scale], parameters=parameters)
	end

	function build_camera_matrix_conic_system(lines_values)
		input_count = size(lines_values)[1]
		rotation_equations_count = 3
		translation_equations_count = 3
		intrinsic_rotation_equations_count = input_count - rotation_equations_count
		@var x y z
		@var f
		@var tx ty tz
		@var lines[1:input_count, 1:3] points_at_infinity[1:input_count, 1:3] dual_quadrics[1:translation_equations_count, 1:4, 1:4]
		intrinsic = build_intrinsic_matrix(f)
		Rₚ = build_rotation_matrix(x, y, z)
		intrinsic_topleft = intrinsic[1:3, 1:3]
		P = build_camera_matrix(intrinsic, Rₚ, [tx; ty; tz]; use_rotation_as_is = true)

		system_to_solve = []
		for i in 1:input_count
			equation = lines[i, :]' * intrinsic_topleft * Rₚ * points_at_infinity[i, :]
			push!(system_to_solve, equation)
		end
		for i in 1:translation_equations_count
			equation = lines[i, :]' * P * dual_quadrics[i, :, :] * P' * lines[i, :]
			push!(system_to_solve, equation)
		end
		parameters::Vector{HomotopyContinuation.ModelKit.Variable} = stack_homotopy_parameters(lines, points_at_infinity, dual_quadrics)
		return System(system_to_solve, variables=[x, y, z, f, tx, ty, tz], parameters=parameters)
	end
end