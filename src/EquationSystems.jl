module EquationSystems
	export stack_homotopy_parameters, build_intrinsic_rotation_conic_system, build_intrinsic_rotation_translation_conic_system, build_camera_matrix_conic_system

	using ..Camera: IntrinsicParameters, build_intrinsic_matrix, build_camera_matrix
	using ..Space: build_rotation_matrix

	using HomotopyContinuation
	using LinearAlgebra: det, I

	module Problems
		using ....Camera: CameraProperties
		mutable struct CylinderCameraContoursProblem
			camera::CameraProperties
			lines::Array{Float64, 2}
			points_at_infinity::Array{Float64, 2}
			dualquadrics::Array{Float64, 3}
		end
	end

	function stack_homotopy_parameters(parameters...)
		stacked_parameters = []
		for parameter in parameters
			stacked_parameters = vcat(stacked_parameters, vec(parameter))
		end
		return stacked_parameters
	end

	function add_rotation_constraints!(system_to_solve, R)
		push!(system_to_solve, det(R) - 1)
		for eq in vec(R * R' - I)
			push!(system_to_solve, eq)
		end
	end

	function build_intrinsic_rotation_conic_system(problems::Vector{Problems.CylinderCameraContoursProblem})
		problems_count = length(problems)
		if problems_count == 0
			throw(ArgumentError("At least one problem is needed"))
		end

		@var scale fₓ fᵧ skew cₓ cᵧ

		system_to_solve = [scale * fₓ * fᵧ - 1]
		variables::Vector{HomotopyContinuation.ModelKit.Variable} = [scale, fₓ, fᵧ, skew, cₓ, cᵧ]
		parameters::Vector{HomotopyContinuation.ModelKit.Variable} = []

		for (index, problem) in enumerate(problems)
			lines_count = size(problem.lines)[1]
			Rparams = [
				Variable("R$(index)", i) for i in 1:3
			]
			R = build_rotation_matrix(Rparams..., false)
			lines = reshape([
				Variable("lines$(index)", i, j)
				for i in 1:lines_count, j in 1:3
			], lines_count, 3)
			points_at_infinity = reshape([
				Variable("points_at_infinity$(index)", i, j)
				for i in 1:lines_count, j in 1:3
			], lines_count, 3)
			intrinsic_topleft = build_intrinsic_matrix(IntrinsicParameters(
				focal_length_x = fₓ,
				focal_length_y = fᵧ,
				principal_point_x = cₓ,
				principal_point_y = cᵧ,
				skew = skew,
			))[1:3, 1:3]

			for line_index in 1:lines_count
				equation = scale * lines[line_index, :]' * intrinsic_topleft * R * points_at_infinity[line_index, :]
				push!(system_to_solve, equation)
			end

			variables = stack_homotopy_parameters(variables, Rparams)
			parameters = stack_homotopy_parameters(parameters, lines, points_at_infinity)
		end

		return System(system_to_solve, variables = variables, parameters = parameters)
	end

	function build_intrinsic_rotation_translation_conic_system(problem::Problems.CylinderCameraContoursProblem)
		lines_count = 3
		@var scale tx ty tz
		@var lines[1:lines_count, 1:3] dual_quadrics[1:lines_count, 1:4, 1:4]
		P = build_camera_matrix(
			problem.camera.intrinsic, 
			problem.camera.quaternion_rotation,
			[tx; ty; tz]
		)

		system_to_solve = [scale - 1.0]
		for i in 1:lines_count
				equation = lines[i, :]' * scale * P * dual_quadrics[i, :, :] * P' * lines[i, :]
				push!(system_to_solve, equation)
		end
		parameters::Vector{HomotopyContinuation.ModelKit.Variable} = stack_homotopy_parameters(lines, dual_quadrics)
		return System(system_to_solve, variables=[scale, tx, ty, tz], parameters=parameters)
	end

	function build_camera_matrix_conic_system(lines_values)
		input_count = size(lines_values)[1]
		rotation_equations_count = 3
		translation_equations_count = 3
		intrinsic_rotation_equations_count = input_count - rotation_equations_count
		@var R[1:3, 1:3]
		@var f
		@var t[1:3]
		@var lines[1:input_count, 1:3] points_at_infinity[1:input_count, 1:3] dual_quadrics[1:translation_equations_count, 1:4, 1:4]
		intrinsic = build_intrinsic_matrix(f)
		P = build_camera_matrix(intrinsic, R, t; use_rotation_as_is = true)

		system_to_solve = []
		add_rotation_constraints!(system_to_solve, R)
		for i in 1:input_count
			equation = lines[i, :]' * intrinsic * R * points_at_infinity[i, :]
			push!(system_to_solve, equation)
		end
		for i in 1:translation_equations_count
			equation = lines[i, :]' * P * dual_quadrics[i, :, :] * P' * lines[i, :]
			push!(system_to_solve, equation)
		end
		variables::Vector{HomotopyContinuation.ModelKit.Variable} = stack_homotopy_parameters([f], R, t)
		parameters::Vector{HomotopyContinuation.ModelKit.Variable} = stack_homotopy_parameters(lines, points_at_infinity, dual_quadrics)
		return System(system_to_solve, variables=variables, parameters=parameters)
	end

	module Minimization
		using ..EquationSystems
		using ..Problems

		using HomotopyContinuation
		using LinearAlgebra: norm

		function bestsolution(system, solutions, parameters)
			best_solution = nothing
			best_residual = Inf
			for solution in solutions
				residual = norm(system(solution, parameters))
				if residual < best_residual
					best_residual = residual
					best_solution = solution
				end
			end
			return best_solution
		end
		function build_intrinsic_rotation_conic_system(problems::Vector{Problems.CylinderCameraContoursProblem})
			sys = EquationSystems.build_intrinsic_rotation_conic_system(problems)
			expressions = sys.expressions
			minimizer = sum(expressions.^2)
			variables = sys.variables
			parameters = sys.parameters

			diff = differentiate(minimizer, variables)

			return System(diff, variables = variables, parameters = parameters)
		end
	end

	module SingleProblem
		function build_intrinsic_rotation_conic_system(lines_values::Matrix{<:Number})
		end
		function build_intrinsic_rotation_translation_conic_system(intrinsic, rotation, lines_values)
		end
	end
end