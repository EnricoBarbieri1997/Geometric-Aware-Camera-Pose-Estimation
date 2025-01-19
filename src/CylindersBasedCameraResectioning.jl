module CylindersBasedCameraResectioning
    include("includes.jl")

	using ..Scene: best_overall_solution!, best_overall_solution_by_steps!, create_scene_instances_and_problems, intrinsic_rotation_system_setup, plot_reconstructed_scene
	using ..EquationSystems.Problems.IntrinsicParameters: Configurations as IntrinsicParametersConfigurations
    using ..Plotting: plot_3dcamera, Plot3dCameraInput
	using ..Printing: print_camera_differences

    using HomotopyContinuation, Serialization

    function main()
        intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ
        scene, problems = create_scene_instances_and_problems(;
            number_of_instances=2,
            number_of_cylinders=2,
            random_seed=14,
            intrinsic_configuration,
        )

        rotation_intrinsic_system, parameters = intrinsic_rotation_system_setup(problems)

        solver = nothing
        starts = nothing
        if true
            solver, starts = solver_startsolutions(
                rotation_intrinsic_system,
                start_system = :total_degree;
                target_parameters = parameters,
            )
            serialize("tmp/total_degree_solver_cache.jld", solver)
            serialize("tmp/total_degree_starts_cache.jld", starts)
        else
            solver = deserialize("tmp/total_degree_solver_cache.jld")
            starts = deserialize("tmp/total_degree_starts_cache.jld")
        end

        chunk_size = 500000
        numberof_start_solutions = length(starts)
        display("Number of start solutions: $numberof_start_solutions. Number of iterations needed: $(ceil(Int, numberof_start_solutions / chunk_size))")
        solution_error = Inf
        for start in Iterators.partition(starts, chunk_size)
            result = nothing
            if true
                result = solve(
                    solver,
                    start;
                )
                serialize("tmp/result_cache.jld", result)
            else
                result = deserialize("tmp/result_cache.jld")
            end
            @info result

            solution_error, _ = best_overall_solution!(
                result,
                scene,
                problems;
                start_error=solution_error,
                intrinsic_configuration,
            )
            if solution_error < Inf
                try
                    serialize("tmp/chunk_test/best_problems_for_error$(round(solution_error; digits=2)).jld", problems)
                catch
                    display("Error saving problems")
                end
            end
            if solution_error < 1e-6
                break
            end
        end

        # for (i, instance) in enumerate(scene.instances)
        #     display("View $i Pre refine")
        #     print_camera_differences(instance.camera, problems[i].camera)
        #     display("--------------------")
        # end

        # refine_best_solution!(scene, problems)
        # plot_reconstructed_scene(scene, problems)

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
    end

    function save_solutions()
        scene, problems = create_scene_instances_and_problems(;
            random_seed=81224,
            number_of_cylinders=4,
            number_of_instances=1,
        )

        rotation_intrinsic_system, parameters = intrinsic_rotation_system_setup(problems)

        result = solve(
            rotation_intrinsic_system,
            target_parameters = parameters,
            start_system = :total_degree,
        )
        @info result

        try
            mkdir("tmp")
        catch
        end

        serialize(
            "tmp/intrinsic_rotation_total_degree_solutions.jld",
            ParametersSolutionsPair(
                parameters,
                solutions(result)
            )
        )

        solution_error, all_possible_solutions = best_intrinsic_rotation_system_solution!(result, scene, problems)

        current_best_result_error = Inf
        reference_translation_result = nothing
        for (i, problem) in enumerate(problems)
            translation_system, parameters = intrinsic_rotation_translation_system_setup(problem)

            result = solve(
                translation_system,
                start_system = :total_degree,
                target_parameters = parameters,
            )
            @info result

            solution_error = best_intrinsic_rotation_translation_system_solution!(result, scene, scene.instances[i], problem)
            if (solution_error < current_best_result_error)
                reference_translation_result = result
            end
            plot_3dcamera(Plot3dCameraInput(
                problem.camera.euler_rotation,
                problem.camera.position,
            ), :green)
        end

        serialize(
            "tmp/translation_total_degree_solutions.jld",
            ParametersSolutionsPair(
                parameters,
                solutions(reference_translation_result)
            )
        )

        for problem in problems
            plot_3dcamera(Plot3dCameraInput(
                problem.camera.euler_rotation,
                problem.camera.position,
            ), :green)
        end

        plot_reconstructed_scene(scene, problems)
    end

    function generate_monodromy_solutions()
        scene, problems = create_scene_instances_and_problems(
            random_seed=777,
            number_of_cylinders=4,
            number_of_instances=1,
        )

        display(scene.figure)

        for (j, problem) in enumerate(problems)
            display("--------------- problem $j ---------------")
            camera = scene.instances[j].camera
            for i in 1:4
                begin #asserts
                    lines = problem.lines
                    points_at_infinity = problem.points_at_infinity
                    dualquadrics = problem.dualquadrics
                    R = camera.quaternion_rotation'
                    cameramatrix = camera.matrix
                    intrinsic = camera.intrinsic
                    display("$(lines[i, :]' * intrinsic * R * points_at_infinity[i, :]), (1) line point")
                    display("$(lines[i, :]' * cameramatrix * dualquadrics[i, :, :] * cameramatrix' * lines[i, :]), (2) line quadric")
                end
            end
            display("------------- end problem $j -------------")
        end

        intrinsic_rotation_system = build_intrinsic_rotation_conic_system(problems)
        fₓ, _, _, skew, fᵧ, _, cₓ, cᵧ, _ = vec(scene.instances[1].camera.intrinsic)
        startingsolution = [fₓ, fᵧ, skew, cₓ, cᵧ]
        parameters = []
        for (i, problem) in enumerate(problems)
            camera = scene.instances[i].camera
            rot = Rotations.params(camera.quaternion_rotation')
            rot = rot / rot[1]
            rot = rot[2:end]
            startingsolution = stack_homotopy_parameters(
                startingsolution,
                rot,
            )
            parameters = stack_homotopy_parameters(
                parameters,
                problem.lines,
                problem.points_at_infinity,
            )
        end
        startingsolution = convert(Vector{Float64}, startingsolution)
        parameters = convert(Vector{Float64}, parameters)
        monodromy_solutions = monodromy_solve(
            intrinsic_rotation_system,
            startingsolution,
            parameters;
            max_loops_no_progress = 60,
        )
        display(monodromy_solutions)

        try
            mkdir("tmp")
        catch
        end

        serialize(
            "tmp/intrinsic_rotation_monodromy_solutions.jld",
            ParametersSolutionsPair(
                parameters,
                solutions(monodromy_solutions)
            )
        )

        problems[1].camera.intrinsic = scene.instances[1].camera.intrinsic
        problems[1].camera.quaternion_rotation = scene.instances[1].camera.quaternion_rotation
        translation_system = build_intrinsic_rotation_translation_conic_system(problems[1])
        startingsolution = scene.instances[1].camera.position
        display("startingsolution: $startingsolution")
        parameters = stack_homotopy_parameters(
            problems[1].lines[1:3, :],
            problems[1].dualquadrics[1:3, :, :]
        )
        startingsolution = convert(Vector{Float64}, startingsolution)
        parameters = convert(Vector{Float64}, parameters)
        monodromy_solutions = monodromy_solve(
            translation_system,
            startingsolution,
            parameters
        )
        display(monodromy_solutions)

        serialize(
            "tmp/translation_monodromy_solutions.jld",
            ParametersSolutionsPair(
                parameters,
                solutions(monodromy_solutions)
            )
        )
    end

    export main, save_solutions, generate_monodromy_solutions
end