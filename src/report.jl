module Report
	using ..Scene: SceneData, best_overall_solution!, best_overall_solution_by_steps!, create_scene_instances_and_problems, intrinsic_rotation_system_setup, plot_scene, plot_reconstructed_scene
	using ..EquationSystems.Problems: CylinderCameraContoursProblem
	using ..EquationSystems.Problems.IntrinsicParameters: Configurations as IntrinsicParametersConfigurations
	using ..Plotting: initfigure, plot_3dcamera, Plot3dCameraInput
	using ..Printing: print_camera_differences, print_error_analysis, create_single_noise_result, save_results_to_json
  	using ..Utils: vector_difference, intrinsic_difference, matrix_difference, rotations_difference, translations_difference
	using ..Geometry: TangentLineNotFound
	
	using CSV, Dates, Glob, HomotopyContinuation, Random, Serialization, Tables
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
	struct ReportConfiguration
		seed::Int
		intrinsic_configuration::IntrinsicParametersConfigurations.T
		noise::Float64
		scene::SceneData
		problems::Vector{CylinderCameraContoursProblem}
	end
	struct ReportResult
		runingtime::Float64
		errors::ErrorsReportData
	end
	struct ReportException
		exception::Exception
	end
	struct ReportData
		configuration::ReportConfiguration
		result::Union{ReportResult, ReportException}
	end

	function is_initialized(report)
		return isa(report, ReportData)
	end

	function is_terminated_successfully(report)
		return is_initialized(report) && isa(report.result, ReportResult)
	end
	
	function multiple_seeds_multiple_configuration(;
		noises = nothing,
		output = false,
	)
		Random.seed!(940)

		number_of_seeds = 50
		configurations = [
			# IntrinsicParametersConfigurations.none,
			# IntrinsicParametersConfigurations.fₓ,
			# IntrinsicParametersConfigurations.fₓ_fᵧ,
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
			(IntrinsicParametersConfigurations.fₓ_fᵧ, [
				(2, 2),
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

		noise_values = if isnothing(noises) collect(0.20:0.05:0.55) else noises end
		
		results = []

		if !isdir("./tmp/reports")
			mkdir("./tmp/reports")
		end
		date_string = Dates.format(now(),"yyyy-mm-dd HH-MM")
		output_path = isa(output, String) ? output : "./tmp/reports/$(date_string)"
		save_in_folder = output === true || isa(output, String)
		if save_in_folder && !isdir(output_path)
			mkdir(output_path)
		end

		scenes = Vector{Union{Tuple{SceneData, Vector{CylinderCameraContoursProblem}}, Nothing}}(undef, number_of_seeds)
		fill!(scenes, nothing)
		seeds = Vector{Int}(undef, number_of_seeds)
		current_noise_results = Vector{Union{ReportData, Exception, Nothing}}(undef, length(noise_values))
		fill!(current_noise_results, nothing)

		for configuration in configurations
			possible_scene_configurations = get(cylinder_views_per_config, configuration, [(2, 1)])
			for scene_configuration in possible_scene_configurations
				number_of_cylinders, number_of_instances = scene_configuration
				for (noise_index, noise) in enumerate(noise_values)
					current_seed_index = 1
					while current_seed_index <= number_of_seeds
						display("Seed index: $current_seed_index. Configuration: $configuration. Scene configuration: $scene_configuration. Noise: $noise")
						start = time()
						report_configuration = nothing
						report_result = nothing
						try
							scene, problems = nothing, nothing
							if isnothing(scenes[current_seed_index])
								seed = rand(Int)
								scenes[current_seed_index] = create_scene_instances_and_problems(;
									number_of_instances,
									number_of_cylinders,
									random_seed = seed,
									intrinsic_configuration = configuration,
									noise,
									plot = false,
								)
								seeds[current_seed_index] = seed
							end
							seed = seeds[current_seed_index]
							scene, problems = scenes[current_seed_index]

							report_configuration = ReportConfiguration(
								seed,
								configuration,
								noise,
								scene,
								problems,
							)

							rotation_intrinsic_system, parameters = intrinsic_rotation_system_setup(problems)

							solver, starts = solver_startsolutions(
								rotation_intrinsic_system,
								start_system = :total_degree;
								target_parameters = parameters,
								endgame_options = EndgameOptions(;
									sing_accuracy = 1e100,
								),
								tracker_options = TrackerOptions(;
									automatic_differentiation = 3,
								)
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

							view_errors = Vector{ViewError}(undef, length(scene.instances))
							for i in 1:length(scene.instances)
								original_camera = scene.instances[i].camera
								calculated_camera = problems[i].camera
								view_errors[i] = ViewError(
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
								)
							end
							current_noise_results[noise_index] = ReportData(
								report_configuration,
								ReportResult(
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
								)
							)

							current_seed_index = current_seed_index + 1
						catch e
							Base.showerror(stdout, e)
							Base.show_backtrace(stdout, catch_backtrace())
							if isa(e, TangentLineNotFound)
								scenes[current_seed_index] = nothing
							else
								if isnothing(report_configuration)
									current_noise_results[noise_index] = e
								else
									current_noise_results[noise_index] = ReportData(report_configuration, ReportException(e))
								end
								current_seed_index = current_seed_index + 1
							end
						end

						display("------------------------------------")
					end
					if length(current_noise_results) > 0
						if save_in_folder
							formatted_noise = "noise-" * replace(string(noise), "." => "-")
							serialize("$(output_path)/$(formatted_noise).jls", current_noise_results)
						else
							append!(results, current_noise_results)
						end
					end
					GC.gc()
				end
			end
		end

		if !save_in_folder
			serialize("./tmp/reports/$(date_string).jls", results)
		end
	end

	function merge_reports(reports_path, output_path)
		reports = []
		files = glob(reports_path)
		for report_path in files
			report = deserialize(report_path)
			append!(reports, report)
		end
		serialize(output_path, reports)
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
			if !is_terminated_successfully(report)
				continue
			end
			height += length(report.result.errors.views)
		end
		data = Matrix{Any}(undef, height, length(header))
		row = 1
		for report in reports
			if !is_terminated_successfully(report)
				continue
			end
			for (j, view) in enumerate(report.result.errors.views)
				data[row, 1] = report.configuration.seed
				data[row, 2] = report.configuration.intrinsic_configuration
				data[row, 3] = length(report.configuration.scene.cylinders)
				data[row, 4] = length(report.configuration.scene.instances)
				data[row, 5] = report.configuration.noise
				data[row, 6] = report.result.runingtime
				data[row, 7] = report.result.errors.intrinsic
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

	function report_error_analysis(report_path, noise_steps; number_of_samples=5, output_path=nothing, save_json=false)
		reports = deserialize(report_path)
		errors_all = Array{Union{Float64, Nothing}, 3}(nothing, 6, length(noise_steps), number_of_samples)
		errors_max = zeros(Float64, 6, length(noise_steps))
		errors_mean = zeros(Float64, 6, length(noise_steps))
		sample_counts = zeros(Int, length(noise_steps))
		errored_count = zeros(Int, length(noise_steps))
		errored = []
		for report in reports
			if !is_initialized(report)
				push!(errored, report)
				continue
			end

			index = findfirst(noise_steps .== report.configuration.noise)

			if !is_terminated_successfully(report)
				errored_count[index] += 1
				push!(errored, report)
				continue
			end
			
			n_views = length(report.result.errors.views)

			mean_cameramatrix_error = reduce((tot, view) -> view.cameramatrix + tot, report.result.errors.views; init=0) / n_views
			mean_rotation_error = reduce((tot, view) -> view.rotation + tot, report.result.errors.views; init=0) / n_views
			mean_translation_error = reduce((tot, view) -> view.translation + tot, report.result.errors.views; init=0) / n_views

			sample_index = sample_counts[index] + errored_count[index] + 1
			errors_all[1:3, index, sample_index] = report.result.errors.intrinsic
			errors_all[4, index, sample_index] = mean_cameramatrix_error
			errors_all[5, index, sample_index] = mean_rotation_error
			errors_all[6, index, sample_index] = mean_translation_error
			errors_mean[1:3, index] += report.result.errors.intrinsic
			errors_mean[4, index] += mean_cameramatrix_error
			errors_mean[5, index] += mean_rotation_error
			errors_mean[6, index] += mean_translation_error
			if (norm(errors_max[1:3, index]) < norm(report.result.errors.intrinsic))
				errors_max[1:3, index] = report.result.errors.intrinsic
			end
			if (errors_max[4, index] < mean_cameramatrix_error)
				errors_max[4, index] = mean_cameramatrix_error
			end
			if (errors_max[5, index] < mean_rotation_error)
				errors_max[5, index] = mean_rotation_error
			end
			if (errors_max[6, index] < mean_translation_error)
				errors_max[6, index] = mean_translation_error
			end
			sample_counts[index] += 1
		end

		errors_mean ./= sample_counts'

		header = []
		for (i, noise) in enumerate(noise_steps)
			push!(header, "$noise ($(sample_counts[i])/$number_of_samples)")
		end

		if !isnothing(output_path)
			mean_output_file = output_path * "mean_errors.csv"
			max_output_file = output_path * "max_errors.csv"
			CSV.write(mean_output_file, Tables.table(errors_mean; header); compact=true)
			CSV.write(max_output_file, Tables.table(errors_max; header); compact=true)
		else
			display("Mean errors")
			print_error_analysis(errors_mean; header)
			display("--------------------")
			display("Max errors")
			print_error_analysis(errors_max; header)
		end

		if save_json
			json_results::Vector{Dict} = []
			for i in eachindex(noise_steps)
				success_rate = sample_counts[i] / (sample_counts[i] + errored_count[i])
				delta_f = filter(x -> x !== nothing, errors_all[1, i, :])
				delta_uv = filter(x -> x !== nothing, errors_all[2, i, :])
				delta_skew = filter(x -> x !== nothing, errors_all[3, i, :])
				delta_r = filter(x -> x !== nothing, errors_all[5, i, :])
				delta_t = filter(x -> x !== nothing, errors_all[6, i, :])
				push!(json_results, create_single_noise_result("ours", noise_steps[i], delta_f, delta_uv, delta_skew, success_rate, delta_r, delta_t))
			end
			for report in errored
				push!(json_results, create_single_noise_result("ours", report.configuration.noise, [], [], [], [], [], []))
			end
			save_results_to_json("assets/methods_compare/ours_results.json", json_results)
		end
	end

	function report_error_configuration_analysis(report_path, configurations; number_of_samples=5, output_path=nothing)
		reports = deserialize(report_path)
		errors_max = zeros(Float64, 6, length(configurations))
		errors_mean = zeros(Float64, 6, length(configurations))
		sample_counts = zeros(Int, length(configurations))
		for report in reports
			if !is_terminated_successfully(report)
				continue
			end
			display(configurations)
			display(report.configuration.intrinsic_configuration)
			index = findfirst(configurations .== UInt8(report.configuration.intrinsic_configuration))
			total_rotation_error = reduce((tot, view) -> view.rotation + tot, report.result.errors.views; init=0)
			total_translation_error = reduce((tot, view) -> view.translation + tot, report.result.errors.views; init=0)
			total_cameramatrix_error = reduce((tot, view) -> view.cameramatrix + tot, report.result.errors.views; init=0)
			errors_mean[1:3, index] += report.result.errors.intrinsic
			errors_mean[4, index] += total_cameramatrix_error
			errors_mean[5, index] += total_rotation_error
			errors_mean[6, index] += total_translation_error
			if (norm(errors_max[1:3, index]) < norm(report.result.errors.intrinsic))
				errors_max[1:3, index] = report.result.errors.intrinsic
			end
			if (errors_max[4, index] < total_cameramatrix_error)
				errors_max[4, index] = total_cameramatrix_error
			end
			if (errors_max[5, index] < total_rotation_error)
				errors_max[5, index] = total_rotation_error
			end
			if (errors_max[6, index] < total_translation_error)
				errors_max[6, index] = total_translation_error
			end
			sample_counts[index] += 1
		end

		errors_mean ./= sample_counts'

		header = []
		for (i, configuration) in enumerate(configurations)
			push!(header, "$configuration ($(sample_counts[i])/$number_of_samples)")
		end

		if !isnothing(output_path)
			mean_output_file = output_path * "mean_errors.csv"
			max_output_file = output_path * "max_errors.csv"
			CSV.write(mean_output_file, Tables.table(errors_mean; header); compact=true)
			CSV.write(max_output_file, Tables.table(errors_max; header); compact=true)
		else
			display("Mean errors")
			print_error_analysis(errors_mean; header)
			display("--------------------")
			display("Max errors")
			print_error_analysis(errors_max; header)
		end
	end

	function explore_report(report_path, seed, intrinsic_configuration, cylinder_view_configuration, noise)
		reports = deserialize(report_path)
		for report in reports
			if (report.configuration.seed == seed &&
				Int(report.configuration.intrinsic_configuration) == intrinsic_configuration &&
				length(report.configuration.scene.cylinders) == cylinder_view_configuration[1] &&
				length(report.configuration.scene.instances) == cylinder_view_configuration[2] &&
				report.configuration.noise == noise
			)
				scene = report.configuration.scene
				problems = report.configuration.problems

				scene.figure = initfigure()
				plot_scene(scene, problems; noise=report.configuration.noise)

				display(scene.figure)

				for (i, instance) in enumerate(scene.instances)
					display("View $i")
					print_camera_differences(instance.camera, problems[i].camera)
					display("--------------------")
				end
		
				for problem in problems
					plot_3dcamera(problem.camera, :green)
				end

				plot_reconstructed_scene(scene, problems)
				display(scene.figure)
				break
			end
		end
	end
end