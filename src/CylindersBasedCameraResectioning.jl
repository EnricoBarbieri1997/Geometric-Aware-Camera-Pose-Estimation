module CylindersBasedCameraResectioning
    const GUI_ENABLED = get(ENV, "GUI_ENABLED", "true") == "true"
    include("includes.jl")

	using ..Scene: ParametersSolutionsPair, best_overall_solution!, best_overall_solution_by_steps!, best_intrinsic_rotation_translation_system_solution!, camera_from_solution, create_scene_instances_and_problems, intrinsic_rotation_system_setup, intrinsic_rotation_translation_system_setup, plot_interactive_scene, plot_reconstructed_scene, split_intrinsic_rotation_parameters
	using ..EquationSystems: stack_homotopy_parameters, variables_jacobian_rank, joint_jacobian_rank
    using ..EquationSystems.Problems.IntrinsicParameters: Configurations as IntrinsicParametersConfigurations
    using ..Plotting
	using ..Printing: print_camera_differences
    using ..Camera: build_camera_matrix
    using ..Homotopies: ParameterHomotopy as MyParameterHomotopy

    using HomotopyContinuation, Observables, Random, Serialization
    using Images, ImageIO, ImageFeatures, Colors

    function main()
        intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ
        scene, problems = create_scene_instances_and_problems(;
            number_of_instances=2,
            number_of_cylinders=2,
            random_seed=14,
            intrinsic_configuration,
            noise=0,
        )

        display(scene.figure)

        rotation_intrinsic_system, parameters = intrinsic_rotation_system_setup(problems)

        start_solutions = nothing
        start_parameters = nothing
        try
            parameters_solutions_pair = deserialize("tmp/intrinsic_rotation_monodromy_solutions.jld")
            start_solutions = parameters_solutions_pair.solutions
            start_parameters = parameters_solutions_pair.start_parameters
        catch e
            @error e
            display("No intrinsic-rotation monodromy")
        end

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

    function from_image(image_path::String)
        # --- Step 1: Load image ---
        img = load(image_path)
        grayimg = Gray.(img)

        # --- Step 2: Detect edges and lines ---
        edges = canny(grayimg, (1.0, 1.0))
        lines = hough_transform_standard(edges; stepsize=1, vote_threshold=100)

        # --- Step 3: Display with clickable interaction ---
        fig = Figure(resolution = (800, 600))
        ax = Plotting.Axis(fig[1, 1];
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
        )
        image!(ax, img)

        # Draw detected lines for reference
        for line in lines
            p1, p2 = line.point1, line.point2
            lines!(ax, [p1[2], p2[2]], [p1[1], p2[1]], color=:red)
        end

        picked_lines = Observable{Vector{Tuple{Point2f, Point2f}}}([])
        temp_points = Observable{Vector{Point2f}}([])

        on(events(ax).mousebutton) do event
            if event.button == Mouse.left && event.action == Mouse.press 
                # mouse_pos = to_world(ax, events(ax).mouseposition[])
                mouse_pos = events(ax).mouseposition[]
                push!(temp_points[], Point2f(mouse_pos))
                temp_points[] = temp_points[]  # trigger update

                # Once 2 points are clicked, store them as a pair
                if length(temp_points[]) == 2
                    push!(picked_lines[], (temp_points[][1], temp_points[][2]))
                    temp_points[] = []
                    println("Selected line $(length(picked_lines[]))")
                end
            end
        end

        # Show selected lines in green
        green_lines_x = lift(picked_lines) do pairs
            x = [0.0, 0.0]
            for pair in pairs
                p1, p2 = pair
                push!(x, [p1[2], p2[2]])
            end
            return x
        end
        green_lines_y = lift(picked_lines) do pairs
            y = [0.0, 0.0]
            for pair in pairs
                p1, p2 = pair
                push!(y, [p1[1], p2[1]])
            end
            return y
        end

        lines!(ax, green_lines_x, green_lines_y; color=:green)

        display(fig)

        # --- Step 4: Wait until user selects 3 line pairs ---
        @info "Please select 6 line pairs (click 2 points for each line)."
        while length(picked_lines[]) < 6
            sleep(0.1)
        end

        @info "Collected all 6 line pairs."

        # --- Step 5: Ask for cylinder parameters ---
        radii = Float64[]
        directions = []

        for i in 1:3
            println("Enter radius for cylinder $i:")
            push!(radii, parse(Float64, readline()))

            println("Enter direction vector for cylinder $i (comma separated, e.g. 1.0,0.0,0.0):")
            dir = split(readline(), ",")
            push!(directions, normalize(parse.(Float64, dir)))
        end

        println("\nSummary:")
        for i in 1:3
            println("Cylinder $i:")
            println("  Line Pair: $(picked_lines[][i])")
            println("  Radius: $(radii[i])")
            println("  Direction: $(directions[i])")
        end
    end

    export explore_path, main, monodromy
end