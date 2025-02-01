module Report
	using ..Scene: SceneData, best_overall_solution!, create_scene_instances_and_problems, intrinsic_rotation_system_setup, plot_scene, plot_reconstructed_scene
	using ..EquationSystems.Problems: CylinderCameraContoursProblem
	using ..EquationSystems.Problems.IntrinsicParameters: Configurations as IntrinsicParametersConfigurations
	using ..Plotting: initfigure, plot_3dcamera, Plot3dCameraInput
	using ..Printing: print_camera_differences, print_error_analysis
  	using ..Utils: vector_difference, intrinsic_difference, matrix_difference, rotations_difference, translations_difference
	
	using CSV, Dates, HomotopyContinuation, Random, Serialization, Tables
	using LinearAlgebra: norm
	
	struct ViewError
		rotation::Float64
		translation::Float64
		cameramatrix::Float64
	end
	struct ErrorsReportData
		intrinsic::Vector{Float64}
		intrinsic_matrix::Float64
		views::Vector{ViewError}
	end
	struct ReportData
		seed::Int
		intrinsic_configuration::IntrinsicParametersConfigurations.T
		noise::Float64
		scene::SceneData
		problems::Vector{CylinderCameraContoursProblem}
		runingtime::Float64
		errors::ErrorsReportData
	end
	
	function multiple_seeds_multiple_configuration(;
		seed_index = nothing,
		noises = nothing,
		save_in_folder = false,
	)
		Random.seed!(938)
		seeds = rand(Int, 5)
		if !isnothing(seed_index)
			seeds = [seeds[seed_index]]
		end

		configurations = [
			# IntrinsicParametersConfigurations.none,
			# IntrinsicParametersConfigurations.fₓ,
			IntrinsicParametersConfigurations.fₓ_fᵧ_cₓ_cᵧ,
			# IntrinsicParametersConfigurations.fₓ_fᵧ_skew_cₓ_cᵧ,
		]

		cylinder_views_per_config = Dict([
			(IntrinsicParametersConfigurations.none, [
				(2, 1),
			]),
			(IntrinsicParametersConfigurations.fₓ, [
				(2, 1),
			]),
			(IntrinsicParametersConfigurations.fₓ_fᵧ_cₓ_cᵧ, [
				# (4, 1),
				(3, 2),
				# (2, 4),
			]),
			(IntrinsicParametersConfigurations.fₓ_fᵧ_skew_cₓ_cᵧ, [
				# (4, 1),
				(3, 2),
				# (2, 5),
			]),
		])

		noise_values = if isnothing(noises) collect(0.0:1.0:10.0) else noises end

		results = []

		for seed in seeds
			for configuration in configurations
				possible_scene_configurations = get(cylinder_views_per_config, configuration, [(2, 1)])
				for scene_configuration in possible_scene_configurations
					for noise in noise_values
						display("Seed: $seed. Configuration: $configuration. Scene configuration: $scene_configuration. Noise: $noise")
						start = time()
						try
							number_of_cylinders, number_of_instances = scene_configuration
							scene, problems = create_scene_instances_and_problems(;
								number_of_instances,
								number_of_cylinders,
								random_seed=seed,
								intrinsic_configuration = configuration,
								noise,
							)

							rotation_intrinsic_system, parameters = intrinsic_rotation_system_setup(problems)

							solver, starts = solver_startsolutions(
								rotation_intrinsic_system,
								start_system = :total_degree;
								target_parameters = parameters,
							)

							chunk_size = 500000
							numberof_start_solutions = length(starts)
							display("Number of start solutions: $numberof_start_solutions. Number of iterations needed: $(ceil(Int, numberof_start_solutions / chunk_size))")
							solution_error = Inf
							for start in Iterators.partition(starts, chunk_size)
								result = nothing
								result = solve(
									solver,
									start;
								)
								@info result

								solution_error, _ = best_overall_solution!(
									result,
									scene,
									problems;
									start_error=solution_error,
									intrinsic_configuration=configuration,
								)
								if solution_error < 1e-6
									break
								end
							end

							view_errors = []
							for i in 1:length(scene.instances)
								original_camera = scene.instances[i].camera
								calculated_camera = problems[i].camera
								push!(view_errors, ViewError(
									rotations_difference(
										original_camera.quaternion_rotation,
										calculated_camera.quaternion_rotation
									),
									translations_difference(
										original_camera.position,
										calculated_camera.position
									),
									matrix_difference(
										original_camera.matrix,
										calculated_camera.matrix
									),
								))
							end
							push!(results, ReportData(
								seed,
								configuration,
								noise,
								scene,
								problems,
								time() - start,
								ErrorsReportData(
									intrinsic_difference(
										problems[1].camera.intrinsic,
										scene.instances[1].camera.intrinsic
									),
									matrix_difference(
										problems[1].camera.intrinsic,
										scene.instances[1].camera.intrinsic
									),
									view_errors,
								),
							))
						catch e
							@error e
							stacktrace(catch_backtrace())
							push!(results, e)
						end

						display("------------------------------------")
					end
				end
			end
		end

		if !isdir("./tmp/reports")
			mkdir("./tmp/reports")
		end
		date_string = Dates.format(now(),"dd-mm-yyyy HH-MM")
		if save_in_folder
			mkdir("./tmp/reports/$(date_string)")
			serialize("./tmp/reports/$date_string)/$(seed_index)-$(join(noises, '_')).jls", results)
		else
			serialize("./tmp/reports/$(date_string).jls", results)
		end
	end

	function report_to_csv(report_path, csv_path)
		reports = deserialize(report_path)
		header = [
			"seed",
			"intrinsic_configuration",
			"number_of_cylinders",
			"number_of_views",
			"noise",
			"running time",
			"intrinsic_error",
			"view",
			"rotation_error",
			"translation_error",
			"cameramatrix_error",
		]
		height = 0
		for report in reports
			height += length(report.errors.views)
		end
		data = Matrix{Any}(undef, height, length(header))
		row = 1
		for report in reports
			for (j, view) in enumerate(report.errors.views)
				data[row, 1] = report.seed
				data[row, 2] = report.intrinsic_configuration
				data[row, 3] = length(report.scene.cylinders)
				data[row, 4] = length(report.scene.instances)
				data[row, 5] = report.noise
				data[row, 6] = report.runingtime
				data[row, 7] = report.errors.intrinsic
				data[row, 8] = j
				data[row, 9] = view.rotation
				data[row, 10] = view.translation
				data[row, 11] = view.cameramatrix

				row += 1
			end
		end

		transform= function (col,val)
			if isa(val, Number) && abs(val) < 1
				return round(val, digits=6)
			end
			return val
		end

		CSV.write(csv_path, Tables.table(data; header); compact=true, transform)
	end

	function report_error_analysis(report_path, noise_steps; output_path=nothing)
		reports = deserialize(report_path)
		errors_max = zeros(Float64, 4, length(noise_steps))
		errors_mean = zeros(Float64, 4, length(noise_steps))
		sample_counts = zeros(Int, length(noise_steps))
		for report in reports
			if !isa(report, ReportData)
				continue
			end
			index = findfirst(noise_steps .== report.noise)
			total_cameramatrix_error = reduce((view, tot) -> view.cameramatrix + tot, 0, report.errors.views)
			errors_mean[1:3, index] += report.errors.intrinsic
			errors_mean[4, index] += total_cameramatrix_error
			if (norm(errors_max[1:3, index]) < norm(report.errors.intrinsic))
				errors_max[1:3, index] = report.errors.intrinsic
			end
			if (errors_max[4, index] < total_cameramatrix_error)
				errors_max[4, index] = total_cameramatrix_error
			end
			sample_counts[index] += 1
		end

		display("Errors tot: $errors_mean")
		display("Errors max: $sample_counts")
		errors_mean ./= sample_counts
		display("Errors mean: $errors_mean")

		if !isnothing(output_path)
			mean_output_file = output_path * "mean_errors.csv"
			max_output_file = output_path * "max_errors.csv"
			CSV.write(mean_output_file, Tables.table(errors_mean; noise_steps); compact=true)
			CSV.write(max_output_file, Tables.table(errors_max; noise_steps); compact=true)
		else
			display("Mean errors")
			print_error_analysis(errors_mean; noise_steps)
			display("--------------------")
			display("Max errors")
			print_error_analysis(errors_max; noise_steps)
		end
	end

	function explore_report(report_path, seed, intrinsic_configuration, cylinder_view_configuration, noise)
		reports = deserialize(report_path)
		for report in reports
			if (report.seed == seed &&
				Int(report.intrinsic_configuration) == intrinsic_configuration &&
				length(report.scene.cylinders) == cylinder_view_configuration[1] &&
				length(report.scene.instances) == cylinder_view_configuration[2] &&
				report.noise == noise
			)
				scene = report.scene
				problems = report.problems

				for (i, instance) in enumerate(scene.instances)
					display("View $i")
					print_camera_differences(instance.camera, problems[i].camera)
					display("--------------------")
				end

				scene.figure = initfigure()
				plot_scene(scene, problems; noise=report.noise)
		
				for problem in problems
					plot_3dcamera(Plot3dCameraInput(
						problem.camera.euler_rotation,
						problem.camera.position,
					), :green)
				end
		
				plot_reconstructed_scene(scene, problems)
				display(scene.figure)
				break
			end
		end
	end
end