module CylindersBasedCameraResectioning
    const GUI_ENABLED = get(ENV, "GUI_ENABLED", "true") == "true"
    const ASSERTS_ENABLED = get(ENV, "ASSERTS_ENABLED", "false") == "true"
    const IMAGE_HEIGHT = 1920
    const IMAGE_WIDTH = 1080
    include("includes.jl")

	using ..Scene: ParametersSolutionsPair, averaged_solution!, best_overall_solution!, best_overall_solution_by_steps!, best_intrinsic_rotation_translation_system_solution!, camera_from_solution, create_scene_instances_and_problems, scene_instances_and_problems_from_files, intrinsic_rotation_system_setup, intrinsic_rotation_translation_system_setup, plot_interactive_scene, plot_reconstructed_scene, split_intrinsic_rotation_parameters
	using ..EquationSystems: stack_homotopy_parameters, variables_jacobian_rank, joint_jacobian_rank
    using ..EquationSystems.Problems.IntrinsicParameters: Configurations as IntrinsicParametersConfigurations
    using ..Plotting
	using ..Printing: print_camera_differences, print_relative_motion_errors
    using ..Camera: build_camera_matrix
    using ..Homotopies: ParameterHomotopy as MyParameterHomotopy

    using HomotopyContinuation, Observables, Random, Serialization

    function main()
        intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ_cₓ_cᵧ
        scene, problems = create_scene_instances_and_problems(;
            number_of_instances=2,
            number_of_cylinders=3,
            random_seed=27,
            intrinsic_configuration,
            noise=0.5,
        )

        display(scene.figure)

        rotation_intrinsic_system, parameters = intrinsic_rotation_system_setup(
            problems;
            algebraic_estimation=true,
        )

        start_solutions = nothing
        start_parameters = nothing
        # try
        #     parameters_solutions_pair = deserialize("tmp/intrinsic_rotation_monodromy_solutions.jld")
        #     start_solutions = parameters_solutions_pair.solutions
        #     start_parameters = parameters_solutions_pair.start_parameters
        # catch e
        #     @error e
        #     display("No intrinsic-rotation monodromy")
        # end

        solver = starts = nothing

        if isnothing(start_parameters)
            solver, starts = solver_startsolutions(
                rotation_intrinsic_system;
                target_parameters = parameters,
                start_system = :total_degree,
            )
        else
            display("Using parameter homotopy")
            geometric_homotopy = MyParameterHomotopy(rotation_intrinsic_system, start_parameters, parameters)
            solver, starts = solver_startsolutions(
                geometric_homotopy,
                start_solutions;
            )
        end

        chunk_size = 500000
        numberof_start_solutions = length(starts)
        display("Number of start solutions: $numberof_start_solutions. Number of iterations needed: $(ceil(Int, numberof_start_solutions / chunk_size))")
        start_error = Inf
        for start in Iterators.partition(starts, chunk_size)
            result = solve(
                solver,
                start;
            )
            # @info result

            start_error, _ = best_overall_solution_by_steps!(
                result,
                problems;
                intrinsic_configuration,
                start_error=start_error,
                scene,
            )
        end

        for (i, instance) in enumerate(scene.instances)
            display("View $i")
            print_camera_differences(instance.camera, problems[i].camera)
            display("--------------------")
        end

        for problem in problems
            plot_3dcamera(Plot3dCameraInput(
                problem.camera.euler_rotation,
                problem.camera.position,
            ), :green)
        end

        plot_reconstructed_scene(scene, problems)

        display(scene.figure)
    end

    function monodromy()
        paramters_solutions_pairs::Vector{ParametersSolutionsPair} = []
        Random.seed!(67)

        for try_index in 1:1000
            try
                intrinsic_rotations_monodromy_solutions::Vector{Vector{ComplexF64}} = []
                translations_monodromy_solutions::Vector{Vector{ComplexF64}} = []

                camera_random_seed = rand(Int)
                scene, problems = create_scene_instances_and_problems(;
                    number_of_instances=2,
                    number_of_cylinders=3,
                    random_seed=camera_random_seed,
                    cylinders_random_seed=14,
                    intrinsic_configuration=IntrinsicParametersConfigurations.fₓ_fᵧ_cₓ_cᵧ,
                )

                rotation_intrinsic_system, parameters = intrinsic_rotation_system_setup(problems)
                # _, total_degree_solutions = total_degree(rotation_intrinsic_system; target_parameters = parameters)

                fₓ, _, _, skew, fᵧ, _, cₓ, cᵧ, _ = vec(scene.instances[1].camera.intrinsic)
                startingsolution = [fₓ, fᵧ, cₓ, cᵧ]

                for instance in scene.instances
                    camera = instance.camera
                    rot = Rotations.params(camera.quaternion_rotation')
                    rot = rot / rot[1]
                    rot = rot[2:end]
                    startingsolution = stack_homotopy_parameters(
                        startingsolution,
                        rot,
                    )
                end

                @info "Jacobians rank"
                @info variables_jacobian_rank(rotation_intrinsic_system, startingsolution, parameters)
                @info joint_jacobian_rank(rotation_intrinsic_system, startingsolution, parameters)

                startingsolution = convert(Vector{Float64}, startingsolution)
                # parameters = convert(Vector{Float64}, parameters)

                push!(intrinsic_rotations_monodromy_solutions, startingsolution)

                while true
                    old_solution_count = length(intrinsic_rotations_monodromy_solutions)
                    monodromy_solutions = monodromy_solve(
                        rotation_intrinsic_system,
                        intrinsic_rotations_monodromy_solutions,
                        parameters;
                    )
                    @info monodromy_solutions
                    push!(intrinsic_rotations_monodromy_solutions, solutions(monodromy_solutions)...)
                    intrinsic_rotations_monodromy_solutions = unique(intrinsic_rotations_monodromy_solutions)
                    new_solution_count = length(intrinsic_rotations_monodromy_solutions)

                    try
                        mkdir("tmp")
                    catch
                    end

                    # for (i, problem) in enumerate(problems)
                    #     problem.camera.intrinsic = scene.instances[i].camera.intrinsic
                    #     problem.camera.quaternion_rotation = scene.instances[i].camera.quaternion_rotation
                    #     translation_system, parameters = intrinsic_rotation_translation_system_setup(problem)
                    #     startingsolution = scene.instances[i].camera.position
                    #     startingsolution = convert(Vector{Float64}, startingsolution)
                    #     monodromy_solutions = monodromy_solve(
                    #         translation_system,
                    #         startingsolution,
                    #         parameters
                    #     )
                    #     @info monodromy_solutions
                    #     push!(translations_monodromy_solutions, solutions(monodromy_solutions)...)
                    #     translations_monodromy_solutions = unique(translations_monodromy_solutions)
                    # end

                    if new_solution_count == old_solution_count
                        break
                    end
                end

                for sol in intrinsic_rotations_monodromy_solutions
                    e = evaluate(rotation_intrinsic_system, sol, parameters)
                    tot_error = sum(abs.(e))
                    if tot_error > 1e-6
                        @error "Error: $tot_error"
                    end
                end

                # append!(intrinsic_rotations_monodromy_solutions, collect(total_degree_solutions))

                push!(paramters_solutions_pairs, ParametersSolutionsPair(
                    parameters,
                    intrinsic_rotations_monodromy_solutions
                ))
            catch e
                @error e
                display("No intrinsic-rotation monodromy")
            end
        end

        serialize(
            "tmp/intrinsic_rotation_monodromy_solutions.jld",
            paramters_solutions_pairs
        )

        # serialize(
        #     "tmp/translation_monodromy_solutions.jld",
        #     ParametersSolutionsPair(
        #         parameters,
        #         translations_monodromy_solutions
        #     )
        # )
    end

    function explore_path()
        intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ
        scene_target, problems_target = create_scene_instances_and_problems(;
            number_of_instances=2,
            number_of_cylinders=2,
            random_seed=14,
            intrinsic_configuration,
            plot=false,
        )
        scene_start, problems_start = create_scene_instances_and_problems(;
            number_of_instances=2,
            number_of_cylinders=2,
            random_seed=67,
            cylinders_random_seed=14,
            intrinsic_configuration,
            plot=false,
        )

        rotation_intrinsic_system_start, parameters_start = intrinsic_rotation_system_setup(problems_start)
        rotation_intrinsic_system_target, parameters_target = intrinsic_rotation_system_setup(problems_target)
        homotopy = ParameterHomotopy(rotation_intrinsic_system_start, parameters_start, parameters_target)
        tracker = Tracker(homotopy)

        fₓ, _, _, skew, fᵧ, _, cₓ, cᵧ, _ = vec(scene_start.instances[1].camera.intrinsic)
        startingsolution = [fₓ, fᵧ]
        for instance in scene_start.instances
            camera = instance.camera
            rot = Rotations.params(camera.quaternion_rotation')
            rot = rot / rot[1]
            rot = rot[2:end]
            startingsolution = stack_homotopy_parameters(
                startingsolution,
                rot,
            )
        end

        pts = []
        for (x, t) in iterator(tracker, startingsolution, 1.0, 0.0)
            push!(pts, real.(x))
        end

        current_figure = initfigure()
        slider = add_slider!(;
            start=1,
            stop=length(pts),
            step=1,
        )
        for i in 1:length(scene_start.instances)
            add_camera_rotation_axis!()
        end
        observable_instances = lift(slider.value) do solution_index
            instances = []
            solution = pts[solution_index]
            intrinsic = rotations_solution = intrinsic_correction = nothing
            try
                intrinsic, rotations_solution, intrinsic_correction = split_intrinsic_rotation_parameters(
                    solution,
                    intrinsic_configuration
                )
            catch e
                @error e
            end
            for (i, instance) in enumerate(scene_start.instances)
                push!(instances, deepcopy(instance))
                if solution_index != 1
                    camera = camera_from_solution(
                        intrinsic,
                        rotations_solution,
                        intrinsic_correction,
                        i
                    )
                    problem = deepcopy(problems_start[i])
                    problem.camera = camera
                    # translation_system, parameters = intrinsic_rotation_translation_system_setup(problem)

                    # translation_result = solve(
                    #     translation_system,
                    #     target_parameters = parameters,
                    #     start_system = :total_degree,
                    # )
                    # @info translation_result

                    # best_intrinsic_rotation_translation_system_solution!(
                    #     translation_result,
                    #     scene_start,
                    #     instances[i],
                    #     problem
                    # )
                    instances[i].camera = problem.camera
                end
            end
            instances
        end

        plot_interactive_scene(;
            scene=scene_target,
            problems=problems_start,
            observable_instances,
            figure=current_figure,
        )
        display(current_figure)
    end

    function markers()
        scene, problems = scene_instances_and_problems_from_files(
            "./assets/test_scenes/markers/scene.json",
            "./assets/methods_compare/real/markers.json";
            number_of_instances=2,
        )
        intrinsic_configuration = problems[1].intrinsic_configuration

        rotation_intrinsic_system, parameters = intrinsic_rotation_system_setup(problems)
        # display(parameters)

        display(scene.figure)
        # return

        solver = starts = nothing

        solver, starts = solver_startsolutions(
            rotation_intrinsic_system;
            target_parameters = parameters,
            start_system = :total_degree,
        )

        chunk_size = 500000
        numberof_start_solutions = length(starts)
        display("Number of start solutions: $numberof_start_solutions. Number of iterations needed: $(ceil(Int, numberof_start_solutions / chunk_size))")
        solution_error = Inf
        for start in Iterators.partition(starts, chunk_size)
            result = solve(
                solver,
                start;
            )
            @info result

            solution_error, _ = best_overall_solution_by_steps!(
                result,
                problems;
                start_error=solution_error,
                intrinsic_configuration,
            )
        end

        for (i, instance) in enumerate(scene.instances)
            display("View $i")
            print_camera_differences(instance.camera, problems[i].camera)
            display("--------------------")
        end

        print_relative_motion_errors(scene, problems)

        for problem in problems
            plot_3dcamera(problem.camera, :green)
        end

        plot_reconstructed_scene(scene, problems)

        display(scene.figure)
    end

    export explore_path, main, monodromy, markers
end