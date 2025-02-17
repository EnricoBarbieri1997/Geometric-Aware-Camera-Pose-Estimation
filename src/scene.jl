module Scene

	using ..Geometry: Line, Cylinder as CylinderType, homogeneous_to_line, line_to_homogenous, homogeneous_line_intercept, get_cylinder_contours
	using ..Space: transformation, random_transformation, identity_transformation, build_rotation_matrix
	using ..Camera: CameraProperties, IntrinsicParameters, build_intrinsic_matrix, build_camera_matrix, lookat_rotation
	using ..Printing: print_camera_differences
	using ..Plotting: initfigure, get_or_add_2d_axis!, clean_plots!, plot_2dpoints, plot_line_2d, Plot3dCameraInput, plot_3dcamera, plot_3dcamera_rotation, Plot3dCylindersInput, plot_3dcylinders, plot_2dcylinders
	using ..EquationSystems: stack_homotopy_parameters, build_intrinsic_rotation_conic_system, build_intrinsic_rotation_translation_conic_system, build_camera_matrix_conic_system
	using ..EquationSystems.Problems: CylinderCameraContoursProblem
	using ..EquationSystems.Problems.IntrinsicParameters: Configurations as IntrinsicParametersConfigurations, has as isIntrinsicEnabled
	# using ..EquationSystems.Minimization: build_intrinsic_rotation_conic_system
	using ..Utils
	using ..Scene
	using ..Cylinder: CylinderProperties, standard_and_dual as standard_and_dual_cylinder
	using ..Conic: ConicProperties
	using Serialization
	using LinearAlgebra: cross, diagm, deg2rad, dot, I, normalize, pinv, svd
	using HomotopyContinuation, Observables, Polynomials, Rotations
	using Random
	using GLMakie: Figure
	struct ParametersSolutionsPair
		start_parameters::Vector{Float64}
		solutions::Union{Vector{Vector{ComplexF64}}, Vector{Vector{Vector{ComplexF64}}}}
	end
	@kwdef mutable struct InstanceConfiguration
		camera::CameraProperties = CameraProperties()
		conics::Vector{ConicProperties} = []
		conics_contours::Array{Float64, 3} = Array{Float64}(undef, 0, 0, 0)
	end
	@kwdef mutable struct SceneData
		cylinders::Vector{CylinderProperties} = []
		instances::Vector{InstanceConfiguration} = []
		figure::Figure = initfigure()
	end

	function create_scene_instances_and_problems(;
		random_seed = 7,
		cylinders_random_seed = random_seed,
		number_of_cylinders = 4,
		number_of_instances = 5,
		noise = 0,
		intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ_skew_cₓ_cᵧ,
		plot = true,
	)
			Random.seed!(cylinders_random_seed)

			scene = SceneData()

			cylinders = []
			for i in 1:number_of_cylinders
					cylinder = CylinderProperties()
					position = normalize(rand(Float64, 3)) * rand_in_range(0.0, 10.0)
					rotation = rand_in_range((-90, 90), 3)
					cylinder.euler_rotation = rotation

					cylinder.transform = transformation(position, cylinder.euler_rotation)
					radius = rand_in_range((1, 3), 2)
					cylinder.radiuses = [radius[1], radius[1]] # TODO Support different radiuses for each cylinder

					axis = cylinder.transform * [0; 0; 1; 0]
					axis = axis ./ axis[3]
					axis = axis[1:3]
					cylinder.geometry = CylinderType(
							position,
							cylinder.radiuses[1],
							axis,
					)

					standard, dual, singularpoint = standard_and_dual_cylinder(cylinder.transform, cylinder.radiuses)
					cylinder.matrix = standard
					cylinder.dual_matrix = dual
					cylinder.singular_point = singularpoint

					push!(cylinders, cylinder)

					begin #asserts
							# @assert cylinder.matrix * cylinder.dual_matrix ≃ diagm([1, 1, 0, 1]) "(-1) The dual quadric is correct"

							@assert cylinder.singular_point' * cylinder.matrix * cylinder.singular_point ≃ 0 "(1) Singular point $(1) belongs to the cylinder $(1)"
							dual_singular_plane = inv(cylinder.transform') * reshape([1, 0, 0, -cylinder.radiuses[1]], :, 1)
							@assert (dual_singular_plane' * cylinder.dual_matrix * dual_singular_plane) ≃ 0 "(2) Perpendicular plane $(1) belongs to the dual cylinder $(1)"

							@assert (cylinder.matrix * cylinder.singular_point) ≃ [0, 0, 0, 0] "(6) Singular point is right null space of cylinder matrix $(i)"

							@assert ((dual_singular_plane' * cylinder.singular_point) ≃ 0 && (dual_singular_plane' * cylinder.dual_matrix * dual_singular_plane) ≃ 0) "(7) Singular plane / point and dual quadric constraints $(i)"
							@assert cylinder.singular_point[4] ≃ 0 "(10) Singular point is at infinity $(i)"
					end
			end

			scene.cylinders = cylinders

			Random.seed!(random_seed)

			instances = []

			focal_length_x = 1
			focal_length_y = 1
			skew = 0
			principal_point_x = 0
			principal_point_y = 0

			if (isIntrinsicEnabled.fₓ(intrinsic_configuration))
					focal_length_x = rand_in_range(2500.0, 3000.0)
			end
			if (isIntrinsicEnabled.fᵧ(intrinsic_configuration))
					focal_length_y = rand_in_range(0.8, 1.0) * focal_length_x
			end
			if (isIntrinsicEnabled.skew(intrinsic_configuration))
					skew = rand_in_range(0, 1)
			end
			if (isIntrinsicEnabled.cₓ(intrinsic_configuration))
					principal_point_x = rand_in_range(1280, 1440)
			end
			if (isIntrinsicEnabled.cᵧ(intrinsic_configuration))
					principal_point_y = principal_point_x * (9/16 + rand_in_range(-0.1, 0.1))
			end
			intrinsic = build_intrinsic_matrix(IntrinsicParameters(;
					focal_length_x,
					focal_length_y,
					principal_point_x,
					principal_point_y,
					skew,
			))

			for i in 1:number_of_instances
					instance = InstanceConfiguration()
					position, rotation_matrix = random_camera_lookingat_center()
					quaternion_camera_rotation = QuatRotation(rotation_matrix)
					euler_rotation = rad2deg.(Rotations.params(RotXYZ(rotation_matrix)))
					camera = CameraProperties(
							position = position,
							euler_rotation = euler_rotation,
							quaternion_rotation = quaternion_camera_rotation,
							intrinsic = intrinsic,
					)
					camera.matrix = build_camera_matrix(intrinsic, quaternion_camera_rotation, position)

					instance.camera = camera
					
					conics = []
					for i in 1:number_of_cylinders
							conic = ConicProperties(
									pinv(camera.matrix') * cylinders[i].matrix * pinv(camera.matrix),
									camera.matrix * cylinders[i].singular_point,
									camera.matrix * cylinders[i].dual_matrix * camera.matrix',
							)
							push!(conics, conic)
					end
					instance.conics = conics

					conics_contours = Array{Float64}(undef, number_of_cylinders, 2, 3)
					for i in 1:number_of_cylinders
							lines = get_cylinder_contours(
									cylinders[i].geometry,
									collect(camera.position),
									camera.matrix
							)
							for (j, line) in enumerate(lines)
									line_homogenous = line_to_homogenous(line)
									conics_contours[i, j, :] = line_homogenous

									begin #asserts
											@assert line_homogenous' * conics[i].dual_matrix * line_homogenous ≃ 0 "(3) Line of projected singular plane $(1) belongs to the dual conic $(1)"
											@assert line_homogenous' * camera.matrix * cylinders[i].singular_point ≃ 0 "(8) Line $(j) of conic $(i) passes through the projected singular point"
											@assert line_homogenous' * camera.matrix * cylinders[i].dual_matrix * camera.matrix' * line_homogenous ≃ 0 "(9) Line $(j) of conic $(i) is tangent to the projected cylinder"
									end
							end
					end
					instance.conics_contours = conics_contours

					push!(instances, instance)
			end

			scene.instances = instances

			intrinsicparamters_count = count_ones(UInt8(intrinsic_configuration))
			problems::Vector{CylinderCameraContoursProblem} = []
			numberoflines_tosolvefor_perinstance = 3 + floor(Int, intrinsicparamters_count/number_of_instances)
			number_of_extra_picks = intrinsicparamters_count % number_of_instances
			for (instance_number, instance) in enumerate(instances)
					conics_contours = instance.conics_contours
					noisy_conic_contours = Array{Float64}(undef, size(conics_contours))
					if noise > 0
						for i in 1:size(conics_contours)[1]
							line1 = conics_contours[i, 1, :]
							line2 = conics_contours[i, 2, :]
							intersection = cross(line1, line2)
							intersection = intersection ./ intersection[3]
							noisy_intersection = intersection + normalize([rand(Float64, 2); 1]) * noise
							noisy_intersection = noisy_intersection ./ noisy_intersection[3]
							noisy_line_1 = cross([0, homogeneous_line_intercept(0, line1), 1], noisy_intersection)
							noisy_line_2 = cross([0, homogeneous_line_intercept(0, line2), 1], noisy_intersection)
							
							noisy_conic_contours[i, 1, :] = noisy_line_1 ./ noisy_line_1[3]
							noisy_conic_contours[i, 2, :] = noisy_line_2 ./ noisy_line_2[3]
						end
					end

					numberoflines_tosolvefor = numberoflines_tosolvefor_perinstance + (instance_number <= number_of_extra_picks ? 1 : 0)

					lines = Matrix{Float64}(undef, numberoflines_tosolvefor, 3)
					noise_free_lines = Matrix{Float64}(undef, numberoflines_tosolvefor, 3)
					points_at_infinity = Matrix{Float64}(undef, numberoflines_tosolvefor, 3)
					dualquadrics = Array{Float64}(undef, numberoflines_tosolvefor, 4, 4)
					possible_picks = collect(1:(number_of_cylinders*2))
					for store_index in (1:numberoflines_tosolvefor)
							line_index = store_index # rand(possible_picks)
							possible_picks = filter(x -> x != line_index, possible_picks)
							i = ceil(Int, line_index / 2)
							j = (line_index - 1) % 2 + 1

							line = conics_contours[i, j, :]
							noise_free_lines[store_index, :] = normalize(line)
							noisy_line = noisy_conic_contours[i, j, :]
							if noise > 0
								lines[store_index, :] = normalize(noisy_line)
							else
								lines[store_index, :] = normalize(line)
							end
							points_at_infinity[store_index, :] = normalize(cylinders[i].singular_point[1:3])
							dualquadrics[store_index, :, :] = cylinders[i].dual_matrix ./ cylinders[i].dual_matrix[4, 4]
					end

					problem = CylinderCameraContoursProblem(
							CameraProperties(),
							lines,
							noise_free_lines,
							points_at_infinity,
							dualquadrics,
							UInt8(intrinsic_configuration),
					)
					push!(problems, problem)
			end

			if plot
				plot_scene(scene, problems; noise)
			end

			return scene, problems
	end

	function plot_scene(scene, problems; noise = 0)
		number_of_cylinders = size(scene.cylinders)[1]
		plot_3dcylinders(Plot3dCylindersInput(
			[cylinder.transform for cylinder in scene.cylinders],
			[cylinder.radiuses for cylinder in scene.cylinders],
			number_of_cylinders,
		))

		for (i, instance) in enumerate(scene.instances)
			camera = instance.camera
			conics = instance.conics
			conics_contours = instance.conics_contours
			plot_3dcamera(Plot3dCameraInput(
					camera.euler_rotation,
					camera.position,
			))
			get_or_add_2d_axis!(i)
			plot_2dpoints([(conic.singular_point ./ conic.singular_point[3])[1:2] for conic in conics]; axindex = i)
			plot_2dcylinders(conics_contours, alpha=0.5; axindex = i)
		end
		if (noise > 0)
			for (i, problem) in enumerate(problems)
				ordered_contours = zeros(size(scene.cylinders)[1] * 2, 3)
				is_first_line = fill(true, size(scene.cylinders)[1])
				for j in 1:size(problem.lines)[1]
					cylinder_index = findfirst((cylinder) -> normalize(cylinder.singular_point[1:3]) == problem.points_at_infinity[j, :], scene.cylinders)
					new_line_index = cylinder_index * 2 - (is_first_line[cylinder_index] ? 1 : 0)
					is_first_line[cylinder_index] = false
					ordered_contours[new_line_index, :, :] = problem.lines[j, :, :]
				end
				noisy_contours = vcat(ordered_contours)
				noisy_contours = reshape(noisy_contours, 2, number_of_cylinders, 3)
				noisy_contours = permutedims(noisy_contours, (2,1,3))
				plot_2dcylinders(noisy_contours; linestyle=:dashdotdot, axindex = i)
			end
		end
	end

	function plot_interactive_scene(;
			scene,
			problems,
			noise = 0,
			observable_instances,
			figure,
		)
		on(observable_instances) do instances
			try
				observable_scene = SceneData(;
					figure,
					cylinders = scene.cylinders,
					instances,
				)
				clean_plots!()
				plot_scene(observable_scene, problems; noise)
				for (i, instance) in enumerate(instances)
					plot_3dcamera_rotation(Plot3dCameraInput(
							scene.instances[i].camera.euler_rotation,
							scene.instances[i].camera.position,
					); axindex = i)
					plot_3dcamera_rotation(Plot3dCameraInput(
							instance.camera.euler_rotation,
							instance.camera.position,
					); color=:green, axindex = i)
				end
			catch e
				@error e
			end
		end
	end

	function split_intrinsic_rotation_parameters(
			solution,
			intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ_skew_cₓ_cᵧ
	)
			intrinsicparamters_count = count_ones(UInt8(intrinsic_configuration))

			intrinsics_solution = solution[1:intrinsicparamters_count]
			focal_length_x = focal_length_y = 1
			principal_point_x = principal_point_y = skew = 0
			intrinsic_solution_index = 1
			if (isIntrinsicEnabled.fₓ(intrinsic_configuration))
					focal_length_x = intrinsics_solution[intrinsic_solution_index]
					intrinsic_solution_index += 1
			end
			if (isIntrinsicEnabled.fᵧ(intrinsic_configuration))
					focal_length_y = intrinsics_solution[intrinsic_solution_index]
					intrinsic_solution_index += 1
			end
			if (isIntrinsicEnabled.skew(intrinsic_configuration))
					skew = intrinsics_solution[intrinsic_solution_index]
					intrinsic_solution_index += 1
			end
			if (isIntrinsicEnabled.cₓ(intrinsic_configuration))
					principal_point_x = intrinsics_solution[intrinsic_solution_index]
					intrinsic_solution_index += 1
			end
			if (isIntrinsicEnabled.cᵧ(intrinsic_configuration))
					principal_point_y = intrinsics_solution[intrinsic_solution_index]
					intrinsic_solution_index += 1
			end

			# Spurious solutions
			if (focal_length_x ≃ 0 || focal_length_y ≃ 0)
				throw(ArgumentError("Spurious solution"))
			end

			intrinsic = build_intrinsic_matrix(IntrinsicParameters(
					focal_length_x = focal_length_x,
					focal_length_y = focal_length_y,
					principal_point_x = principal_point_x,
					principal_point_y = principal_point_y,
					skew = skew,
			))
			intrinsic_correction = I
			if (focal_length_x < 0)
					intrinsic_correction *= [
							-1 0 0;
							0 1 0;
							0 0 1;
					]
			end
			if (focal_length_y < 0 && skew <= 0)
					intrinsic_correction *= [
							1 0 0;
							1 -1 0;
							0 0 1;
					]
			end
			intrinsic_solution = intrinsic * intrinsic_correction
			rotations_solution = solution[(intrinsicparamters_count + 1):end]

			return intrinsic_solution, rotations_solution, intrinsic_correction
	end

	function camera_from_solution(
		intrinsic,
		rotations_solution,
		intrinsic_correction,
		index,
	)
		quat = [1; rotations_solution[(index-1)*3+1:index*3]]
		quat = quat / norm(quat)
		camera_extrinsic_rotation = QuatRotation(quat) * inv(intrinsic_correction)

		return CameraProperties(
				euler_rotation = rad2deg.(eulerangles_from_rotationmatrix(camera_extrinsic_rotation')),
				quaternion_rotation = camera_extrinsic_rotation',
				intrinsic = intrinsic,
		)
	end

	function intrinsic_rotation_system_setup(problems)
			rotation_intrinsic_system = build_intrinsic_rotation_conic_system(
				problems;
				focal_guess = 2000.0
			)
			parameters = []
			for problem in problems
				parameters = stack_homotopy_parameters(
					parameters,
					problem.lines,
				)
			end
			parameters = convert(Vector{Float64}, parameters)

			return rotation_intrinsic_system, parameters
	end

	function best_intrinsic_rotation_system_solution!(
			result,
			scene,
			problems;
			start_error = Inf,
			intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ_skew_cₓ_cᵧ
	)
			solution_error = start_error
			solutions_to_try = real_solutions(result)
			all_possible_solutions = []
			best_solution = nothing
			for solution in solutions_to_try
					intrinsic = rotations_solution = intrinsic_correction = nothing
					try
						intrinsic, rotations_solution, intrinsic_correction = split_intrinsic_rotation_parameters(
							solution,
							intrinsic_configuration
						)
					catch
						continue
					end

					current_error = 0
					possible_cameras = []
					individual_problem_max_error = -Inf
					for i in 1:length(problems)
							individual_problem_error = 0
							quat = [1; rotations_solution[(i-1)*3+1:i*3]]
							quat = quat / norm(quat)
							camera_extrinsic_rotation = QuatRotation(quat) * inv(intrinsic_correction)

							possible_camera = camera_from_solution(
								intrinsic,
								rotations_solution,
								intrinsic_correction,
								i
							)
							push!(possible_cameras, possible_camera)

							for (j, contour) in enumerate(eachslice(scene.instances[i].conics_contours, dims=1))
									for line in eachslice(contour, dims=1)
											eq = line' * intrinsic * camera_extrinsic_rotation * scene.cylinders[j].singular_point[1:3]
											individual_problem_error += abs(eq)
									end
							end
							if (individual_problem_error > individual_problem_max_error)
									individual_problem_max_error = individual_problem_error
							end
					end
					current_error = individual_problem_max_error
					push!(all_possible_solutions, possible_cameras[1])

					if (current_error < solution_error)
							solution_error = current_error
							best_solution = solution
							for (i, problem) in enumerate(problems)
									problem.camera = possible_cameras[i]
							end
					end
			end

			if (!isnothing(best_solution))
				paths = path_results(result)
				best_path = paths[findall(x -> real.(x.solution) == best_solution, paths)[1]]
				display("The best starting solution was $(best_path.start_solution)")
			end

			return solution_error, all_possible_solutions
	end

	function intrinsic_rotation_translation_system_setup(problem)
			translation_system = build_intrinsic_rotation_translation_conic_system(
					problem
			)
			parameters = stack_homotopy_parameters(problem.lines[1:3, :])

			return translation_system, parameters
	end
	function best_intrinsic_rotation_translation_system_solution!(
			result,
			scene,
			instance,
			problem,
	)
			solution_error = Inf
			solutions_to_try = real_solutions(result)
			reference_translation_result = nothing
			for solution in solutions_to_try
					tx, ty, tz = solution
					position = [tx, ty, tz]

					camera_matrix = build_camera_matrix(
							problem.camera.intrinsic,
							problem.camera.quaternion_rotation,
							position
					)

					acceptable = true
					current_error = 0
					for (i, contour) in enumerate(eachslice(instance.conics_contours, dims=1))
							for line in eachslice(contour, dims=1)
									eq = line' * camera_matrix * scene.cylinders[i].dual_matrix * camera_matrix' * line
									current_error += abs(eq)
							end
					end
					if (current_error < solution_error)
							solution_error = current_error
							problem.camera.position = position
							problem.camera.matrix = camera_matrix
					end
			end

			return solution_error
	end

	function best_overall_solution!(
			result,
			scene,
			problems;
			start_error = Inf,
			intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ_skew_cₓ_cᵧ
	)
			solution_error = start_error
			solutions_to_try = real_solutions(result)
			all_possible_solutions = []
			for solution in solutions_to_try
					intrinsic = rotations_solution = intrinsic_correction = nothing
					try
						intrinsic, rotations_solution, intrinsic_correction = split_intrinsic_rotation_parameters(
							solution,
							intrinsic_configuration
						)
					catch
						continue
					end

					current_error = 0
					possible_cameras = []
					# individual_problem_min_error = Inf
					individual_problem_max_error = -Inf
					for (i, problem) in enumerate(problems)
							individual_problem_error = 0
							quat = [1; rotations_solution[(i-1)*3+1:i*3]]
							quat = quat / norm(quat)
							camera_extrinsic_rotation = QuatRotation(quat) * inv(intrinsic_correction)

							problem_upto_translation = CylinderCameraContoursProblem(
									CameraProperties(
											euler_rotation = rad2deg.(eulerangles_from_rotationmatrix(camera_extrinsic_rotation')),
											quaternion_rotation = camera_extrinsic_rotation',
											intrinsic = intrinsic,
									),
									problem.lines,
									problem.noise_free_lines,
									problem.points_at_infinity,
									problem.dualquadrics,
									problem.intrinsic_configuration,
							)

							translation_system, parameters = intrinsic_rotation_translation_system_setup(problem_upto_translation)

							start_solutions = nothing
							start_parameters = nothing
							try
									parameters_solutions_pair = deserialize("tmp/translation_monodromy_solutions.jld")
									start_solutions = parameters_solutions_pair.solutions
									start_parameters = parameters_solutions_pair.start_parameters
							catch
									display("No translation monodromy")
							end

							translation_result = solve(
									translation_system,
									start_solutions;
            			start_parameters,
									target_parameters = parameters,
									start_system = :total_degree,
							)
							@info result

							individual_problem_error += best_intrinsic_rotation_translation_system_solution!(
									translation_result,
									scene,
									scene.instances[i],
									problem_upto_translation
							)

							possible_camera = problem_upto_translation.camera
							push!(possible_cameras, possible_camera)

							for (i, contour) in enumerate(eachslice(scene.instances[i].conics_contours, dims=1))
									for line in eachslice(contour, dims=1)
											eq = line' * intrinsic * camera_extrinsic_rotation * scene.cylinders[i].singular_point[1:3]
											individual_problem_error += abs(eq)
									end
							end
							# if (individual_problem_error < individual_problem_min_error)
							#     individual_problem_min_error = individual_problem_error
							# end
							if (individual_problem_error > individual_problem_max_error)
									individual_problem_max_error = individual_problem_error
							end
					end
					current_error = individual_problem_max_error
					push!(all_possible_solutions, possible_cameras[1])

					if (current_error < solution_error)
							solution_error = current_error
							for (i, problem) in enumerate(problems)
									problem.camera = possible_cameras[i]
							end
					end
			end

			return solution_error, all_possible_solutions
	end

	function best_overall_solution_by_steps!(
			result,
			scene,
			problems;
			start_error = Inf,
			intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ_skew_cₓ_cᵧ
	)
			solution_error, all_possible_solutions = best_intrinsic_rotation_system_solution!(
					result,
					scene,
					problems;
					start_error=start_error,
					intrinsic_configuration,
			)

			for (i, problem) in enumerate(problems)
					display("Problem $i")
					translation_system, parameters = intrinsic_rotation_translation_system_setup(problem)

					start_solutions = nothing
					start_parameters = nothing
					try
							parameters_solutions_pair = deserialize("tmp/translation_monodromy_solutions.jld")
							start_solutions = parameters_solutions_pair.solutions
							start_parameters = parameters_solutions_pair.start_parameters
					catch e
							display(e)
							display("No translation monodromy")
					end

					result = solve(
							translation_system,
							start_solutions;
							start_parameters,
							target_parameters = parameters,
							start_system = :total_degree,
					)
					@info result

					solution_error += best_intrinsic_rotation_translation_system_solution!(result, scene, scene.instances[i], problem)
					plot_3dcamera(Plot3dCameraInput(
							problem.camera.euler_rotation,
							problem.camera.position,
					), :green)
			end

			return solution_error, all_possible_solutions
	end

	function refine_best_solution!(
			scene,
			problems
	)
			display("Refine step")
			# intrinsic_inverse = inv(problems[1].camera.intrinsic)
			refine_problems = [
					CylinderCameraContoursProblem(
							CameraProperties(
									intrinsic = problems[1].camera.intrinsic,
							),
							problem.lines,
							problem.noise_free_lines,
							# vcat([(intrinsic_inverse * line)' for line in eachrow(problem.lines)]...),
							# vcat([(intrinsic_inverse * line)' for line in eachrow(problem.noise_free_lines)]...),
							problem.points_at_infinity,
							problem.dualquadrics,
							UInt8(IntrinsicParametersConfigurations.none),
					)
					for problem in problems
			]
			intrinsic_rotation_system, parameters = intrinsic_rotation_system_setup(refine_problems)
			# start_solutions = stack_homotopy_parameters(
			#     [Rotations.params(problem.camera.quaternion_rotation)[2:4] for problem in problems]...
			# )
			@info intrinsic_rotation_system
			result = solve(
					intrinsic_rotation_system,
					# start_solutions;
					target_parameters = parameters,
					start_system = :total_degree,
			)
			display("---------------------------")
			@info result
			best_overall_solution_by_steps!(
					result,
					scene,
					refine_problems;
					intrinsic_configuration=UInt8(IntrinsicParametersConfigurations.none),
			)
			for (i, problem) in enumerate(problems)
					problem.camera.euler_rotation = refine_problems[i].camera.euler_rotation
					problem.camera.quaternion_rotation = refine_problems[i].camera.quaternion_rotation
					problem.camera.position = refine_problems[i].camera.position
					problem.camera.matrix = convert(Matrix{Float64}, build_camera_matrix(
							problem.camera.intrinsic,
							problem.camera.quaternion_rotation,
							problem.camera.position
					))
			end
	end

	function plot_reconstructed_scene(scene, problems)
			number_of_cylinders = size(scene.cylinders)[1]
			for (i, problem) in enumerate(problems)
					reconstructued_contours = Array{Float64}(undef, number_of_cylinders, 2, 3)
					for i in 1:number_of_cylinders
							lines = get_cylinder_contours(
									scene.cylinders[i].geometry,
									collect(problem.camera.position),
									problem.camera.matrix
							)
							for (j, line) in enumerate(lines)
									line_homogenous = line_to_homogenous(line)
									reconstructued_contours[i, j, :] = line_homogenous
							end
					end

					plot_2dcylinders(reconstructued_contours, linestyle=:dash; axindex = i)
			end
	end
end