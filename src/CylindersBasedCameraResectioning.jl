module CylindersBasedCameraResectioning
    include("includes.jl")

	using ..Scene: ParametersSolutionsPair, best_overall_solution!, best_overall_solution_by_steps!, create_scene_instances_and_problems, intrinsic_rotation_system_setup, intrinsic_rotation_translation_system_setup, plot_reconstructed_scene
	using ..EquationSystems: stack_homotopy_parameters
    using ..EquationSystems.Problems.IntrinsicParameters: Configurations as IntrinsicParametersConfigurations
    using ..Plotting: plot_3dcamera, Plot3dCameraInput
	using ..Printing: print_camera_differences

    using HomotopyContinuation, Random, Serialization

    function main()
        intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ
        scene, problems = create_scene_instances_and_problems(;
            number_of_instances=2,
            number_of_cylinders=2,
            random_seed=14,
            intrinsic_configuration,
            noise=0.5,
        )

        display(scene.figure)

        rotation_intrinsic_system, parameters = intrinsic_rotation_system_setup(problems)

        start_solutions = nothing
        start_parameters = nothing
        try
            parameters_solutions_pair = deserialize("tmp/intrinsic_rotation_monodromy_solutions.jld")
            start_solutions = parameters_solutions_pair.solutions
            start_parameters = parameters_solutions_pair.start_parameters
        catch
            display("No intrinsic-rotation monodromy")
        end

        solver, starts = solver_startsolutions(
            rotation_intrinsic_system,
            start_solutions;
            start_parameters,
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
                scene,
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
        intrinsic_rotations_monodromy_solutions::Vector{Vector{ComplexF64}} = []
        translations_monodromy_solutions::Vector{Vector{ComplexF64}} = []

        scene, problems = create_scene_instances_and_problems(
            number_of_instances=2,
            number_of_cylinders=2,
            random_seed=67,
            cylinders_random_seed=14,
            intrinsic_configuration=IntrinsicParametersConfigurations.fₓ_fᵧ,
        )

        rotation_intrinsic_system, parameters = intrinsic_rotation_system_setup(problems)

        fₓ, _, _, skew, fᵧ, _, cₓ, cᵧ, _ = vec(scene.instances[1].camera.intrinsic)
        startingsolution = [fₓ, fᵧ]
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

        startingsolution = convert(Vector{Float64}, startingsolution)
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

            display("$new_solution_count, $old_solution_count")

            if new_solution_count == old_solution_count
                break
            end
        end 

        display(intrinsic_rotations_monodromy_solutions)
        display(translations_monodromy_solutions)
        serialize(
            "tmp/intrinsic_rotation_monodromy_solutions.jld",
            ParametersSolutionsPair(
                parameters,
                intrinsic_rotations_monodromy_solutions
            )
        )

        # serialize(
        #     "tmp/translation_monodromy_solutions.jld",
        #     ParametersSolutionsPair(
        #         parameters,
        #         translations_monodromy_solutions
        #     )
        # )
    end

    export main, monodromy
end