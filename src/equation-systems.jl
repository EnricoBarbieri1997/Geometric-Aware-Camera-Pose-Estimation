module EquationSystems
	export stack_homotopy_parameters, build_intrinsic_rotation_conic_system, build_intrinsic_rotation_translation_conic_system, build_camera_matrix_conic_system

	using ..Camera: IntrinsicParameters, build_intrinsic_matrix, build_camera_matrix
	using ..Space: build_rotation_matrix

	using HomotopyContinuation
	using LinearAlgebra: det, I

	module Problems
		using ....Camera: CameraProperties
		module IntrinsicParameters
			@enum T begin
				focal_length_x = 0b00001
				focal_length_y = 0b00010
				skew = 0b00100
				principal_point_x = 0b01000
				principal_point_y = 0b10000
			end
			function Base.:|(a::T, b::T)
				UInt8(a) | UInt8(b)
			end
			function Base.:&(a::T, b::T)
				UInt8(a) & UInt8(b)
			end
			function Base.:|(a::Any, b::T)
				a | UInt8(b)
			end
			function Base.:&(a::Any, b::T)
				a & UInt8(b)
			end
			function Base.:|(a::T, b::Any)
				UInt8(a) | b
			end
			function Base.:&(a::T, b::Any)
				UInt8(a) & b
			end
			module has
				using ..IntrinsicParameters: focal_length_x, focal_length_y, skew as skewParameter, principal_point_x, principal_point_y
				function fₓ(config::Any)
					return (UInt8(config) & focal_length_x) != 0
				end
				function fᵧ(config::Any)
					return (UInt8(config) & focal_length_y) != 0
				end
				function skew(config::Any)
					return (UInt8(config) & skewParameter) != 0
				end
				function cₓ(config::Any)
					return (UInt8(config) & principal_point_x) != 0
				end
				function cᵧ(config::Any)
					return (UInt8(config) & principal_point_y) != 0
				end
			end
			module Configurations
				using ..IntrinsicParameters: focal_length_x, focal_length_y, skew, principal_point_x, principal_point_y
				@enum T begin
					none = 0
					fₓ = UInt8(focal_length_x)
					fᵧ = UInt8(focal_length_y)
					fₓ_fᵧ = focal_length_x | focal_length_y
					fₓ_fᵧ_skew = focal_length_x | focal_length_y | skew
					fₓ_fᵧ_skew_cₓ = focal_length_x | focal_length_y | skew | principal_point_x
					fₓ_fᵧ_skew_cᵧ = focal_length_x | focal_length_y | skew | principal_point_y
					fₓ_fᵧ_cₓ_cᵧ = focal_length_x | focal_length_y | principal_point_x | principal_point_y
					fₓ_fᵧ_skew_cₓ_cᵧ = focal_length_x | focal_length_y | skew | principal_point_x | principal_point_y
				end
			end
		end
		mutable struct CylinderCameraContoursProblem
			camera::CameraProperties
			lines::Array{Float64, 2}
			noise_free_lines::Array{Float64, 2}
			points_at_infinity::Array{Float64, 2}
			dualquadrics::Array{Float64, 3}
			intrinsic_configuration::UInt8
		end

		function CylinderCameraContoursProblem(
			camera::CameraProperties,
			lines::Array{Float64, 2},
			noise_free_lines::Array{Float64, 2},
			points_at_infinity::Array{Float64, 2},
			dualquadrics::Array{Float64, 3}
		)
			return CylinderCameraContoursProblem(
				camera,
				lines,
				noise_free_lines,
				points_at_infinity,
				dualquadrics,
				IntrinsicParameters.focal_length_x | IntrinsicParameters.focal_length_y | IntrinsicParameters.skew | IntrinsicParameters.principal_point_x | IntrinsicParameters.principal_point_y
			)
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

	function build_intrinsic_rotation_conic_system(
		problems::Vector{Problems.CylinderCameraContoursProblem};
		focal_guess = 1.0	
	)
		problems_count = length(problems)
		if problems_count == 0
			throw(ArgumentError("At least one problem is needed"))
		end

		# @var focal_guess_enforcer

		fₓ = fᵧ = 1
		skew = cₓ = cᵧ = 0

		system_to_solve = []
		variables::Vector{HomotopyContinuation.ModelKit.Variable} = []
		parameters::Vector{HomotopyContinuation.ModelKit.Variable} = []

		intrinsic_configuration = problems[1].intrinsic_configuration
		if Problems.IntrinsicParameters.has.fₓ(intrinsic_configuration)
			@var fₓ
			push!(variables, fₓ)
		end
		if Problems.IntrinsicParameters.has.fᵧ(intrinsic_configuration)
			@var fᵧ
			push!(variables, fᵧ)
		end
		if Problems.IntrinsicParameters.has.skew(intrinsic_configuration)
			@var skew
			push!(variables, skew)
		end
		if Problems.IntrinsicParameters.has.cₓ(intrinsic_configuration)
			@var cₓ
			push!(variables, cₓ)
		end
		if Problems.IntrinsicParameters.has.cᵧ(intrinsic_configuration)
			@var cᵧ
			push!(variables, cᵧ)
		end

		intrinsic = build_intrinsic_matrix(IntrinsicParameters(
			focal_length_x = fₓ^2,
			focal_length_y = fᵧ^2,
			principal_point_x = cₓ,
			principal_point_y = cᵧ,
			skew = skew^2,
		))

		if UInt8(intrinsic_configuration) == 0
			intrinsic = problems[1].camera.intrinsic
		end

		noises = []
		for (index, problem) in enumerate(problems)
			noise = Variable("noise", index)
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

			for line_index in 1:lines_count
				equation = lines[line_index, :]' * intrinsic * R * points_at_infinity[line_index, :]
				push!(system_to_solve, equation)
			end

			variables = stack_homotopy_parameters(variables, Rparams)
			parameters = stack_homotopy_parameters(parameters, lines, points_at_infinity)
			# push!(noises, noise)
			# push!(system_to_solve, noise)
		end

		# display("focal_guess_enforcer: $focal_guess")
		# push!(system_to_solve, fₓ+fᵧ+focal_guess_enforcer - focal_guess)
		# variables = stack_homotopy_parameters(variables, [focal_guess_enforcer])
		# variables = stack_homotopy_parameters(variables, noises)

		return System(system_to_solve, variables = variables, parameters = parameters)
	end

	function build_intrinsic_rotation_translation_conic_system(problem::Problems.CylinderCameraContoursProblem)
		lines_count = 3
		@var tx ty tz
		@var lines[1:lines_count, 1:3] dual_quadrics[1:lines_count, 1:4, 1:4]
		P = build_camera_matrix(
			problem.camera.intrinsic, 
			problem.camera.quaternion_rotation,
			[tx; ty; tz]
		)

		system_to_solve = []
		for i in 1:lines_count
				equation = lines[i, :]' * P * dual_quadrics[i, :, :] * P' * lines[i, :]
				push!(system_to_solve, equation)
		end
		parameters::Vector{HomotopyContinuation.ModelKit.Variable} = stack_homotopy_parameters(lines, dual_quadrics)
		return System(system_to_solve, variables=[tx, ty, tz], parameters=parameters)
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
		function build_intrinsic_rotation_conic_system(
			problems::Vector{Problems.CylinderCameraContoursProblem};
			args...,
		)
			sys = EquationSystems.build_intrinsic_rotation_conic_system(problems; args...)
			expressions = sys.expressions
			minimizer = sum(expressions.^2)
			variables = sys.variables
			parameters = sys.parameters

			diff = expand.(differentiate(minimizer, variables))

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