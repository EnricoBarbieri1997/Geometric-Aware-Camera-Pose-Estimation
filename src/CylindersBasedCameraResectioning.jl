module CylindersBasedCameraResectioning
    include("includes.jl")

    using .Geometry: Line, Cylinder as CylinderType, line_to_homogenous, get_cylinder_contours
    using .Space: transformation, random_transformation, identity_transformation, build_rotation_matrix
    using .Camera: CameraProperties, IntrinsicParameters, build_intrinsic_matrix, build_camera_matrix, lookat_rotation
    using .Plotting: initfigure, add_2d_axis!, plot_2dpoints, plot_line_2d, Plot3dCameraInput, plot_3dcamera, Plot3dCylindersInput, plot_cylinders_contours, plot_3dcylinders, plot_2dcylinders
    using .EquationSystems: stack_homotopy_parameters, build_intrinsic_rotation_conic_system, build_intrinsic_rotation_translation_conic_system, build_camera_matrix_conic_system
    using .EquationSystems.Problems: CylinderCameraContoursProblem
    using .Debug
    using .Utils
    using LinearAlgebra: diagm, deg2rad, dot, normalize, pinv, svd
    using HomotopyContinuation, Polynomials, Rotations
    using Random
    using Serialization
    using GLMakie: Figure

    struct MonodromyParametersSolutionsPair
        start_parameters::Vector{Float64}
        solutions::Vector{Vector{ComplexF64}}
    end
    @kwdef mutable struct InstanceConfiguration
        camera::CameraProperties = CameraProperties()
        conics::Vector{Conic.ConicProperties} = []
        conics_contours::Array{Float64, 3} = Array{Float64}(undef, 0, 0, 0)
    end
    @kwdef mutable struct Scene
        cylinders::Vector{Cylinder.CylinderProperties} = []
        instances::Vector{InstanceConfiguration} = []
        figure::Figure = initfigure()
    end

    function main()
        scene, problems = create_scene_instances_and_problems(number_of_instances=1)
        camera = scene.instances[1].camera

        parameters_solutions_pair = nothing
        try
            parameters_solutions_pair = deserialize("tmp/intrinsic_rotation_monodromy_solutions.jld")
        catch
            error("generate ir monodromy first")
            return 1
        end

        rotation_intrinsic_system = build_intrinsic_rotation_conic_system(problems)
        parameters = []
        for problem in problems
            parameters = stack_homotopy_parameters(
                parameters,
                problem.lines,
                problem.points_at_infinity,
            )
        end
        parameters = convert(Vector{Float64}, parameters)

        result = solve(
            rotation_intrinsic_system,
            # parameters_solutions_pair.solutions,
            # start_parameters = parameters_solutions_pair.start_parameters,
            target_parameters = parameters
        )
        @info result

        solution_error = Inf
        solutions_to_try = real_solutions(result)
        for solution in solutions_to_try
            intrinsic = build_intrinsic_matrix(IntrinsicParameters(
                focal_length_x = solution[2],
                focal_length_y = solution[3],
                principal_point_x = solution[5],
                principal_point_y = solution[6],
                skew = solution[4],
            )) * solution[1]
            rotations_solution = solution[7:end]

            acceptable = true
            current_error = 0
            possible_cameras = []
            for (i, problem) in enumerate(problems)
                camera_extrinsic_rotation = QuatRotation(
                    1,
                    rotations_solution[(i-1)*3+1:i*3]...
                )

                possible_camera = CameraProperties(
                    euler_rotation = rad2deg.(eulerangles_from_rotationmatrix(camera_extrinsic_rotation')),
                    quaternion_rotation = camera_extrinsic_rotation',
                    intrinsic = intrinsic,
                )
                push!(possible_cameras, possible_camera)

                for (i, contour) in enumerate(eachslice(scene.instances[i].conics_contours, dims=1))
                    for line in eachslice(contour, dims=1)
                        eq = line' * intrinsic[1:3, 1:3] * camera_extrinsic_rotation * scene.cylinders[i].singular_point[1:3]
                        current_error += abs(eq)

                        if (!(eq ≃ 0))
                            acceptable = false
                        end
                    end
                    if (!acceptable)
                        break
                    end
                end
            end

            if (current_error < solution_error)
                solution_error = current_error
                for (i, problem) in enumerate(problems)
                    problem.camera = possible_cameras[i]
                end
            end
        end

        camera_calculated = problems[1].camera

        begin #asserts
            @assert isnothing(camera_calculated) == false  "(11) Found rotation satisfies constraints"
        end

        display("O: $(round.(camera.quaternion_rotation, digits=2))")
        display("C: $(round.(camera_calculated.quaternion_rotation, digits=2))")
        display("CI: $(round.(camera_calculated.quaternion_rotation', digits=2))")

        display("Error of the best rotation solution: $(solution_error)")
        display("Best solution for rotation: $(camera_calculated.euler_rotation)")
        display("Actual rotation: $(camera.euler_rotation)")
        display("Difference between rotations: $(rotations_difference(camera_calculated.quaternion_rotation, camera.quaternion_rotation))")

        display("Calculated intrinsic: $(camera_calculated.intrinsic)")
        display("Actual intrinsic: $(camera.intrinsic)")

        parameters_solutions_pair = nothing
        try
            parameters_solutions_pair = deserialize("tmp/translation_monodromy_solutions.jld")
        catch
            error("generate t monodromy first")
            return 1
        end

        for (i, problem) in enumerate(problems)
            translation_system = build_intrinsic_rotation_translation_conic_system(
                problem
            )
            parameters = stack_homotopy_parameters(problem.lines[1:3, :], problem.dualquadrics[1:3, :, :])

            result = solve(
                translation_system,
                # parameters_solutions_pair.solutions,
                # start_parameters = parameters_solutions_pair.start_parameters,
                target_parameters = parameters,
            )
            @info result

            solution_error = Inf
            solutions_to_try = real_solutions(result)
            for solution in solutions_to_try
                scale, tx, ty, tz = solution
                position = [tx, ty, tz] / scale

                display("$(position) $(scale)")

                camera_matrix = build_camera_matrix(
                    problem.camera.intrinsic,
                    problem.camera.quaternion_rotation,
                    position
                )

                acceptable = true
                current_error = 0
                for (i, contour) in enumerate(eachslice(scene.instances[i].conics_contours, dims=1))
                    for line in eachslice(contour, dims=1)
                        eq = line' * camera_matrix * scene.cylinders[i].dual_matrix * camera_matrix' * line
                        current_error += abs(eq)

                        if (!(eq ≃ 0))
                            acceptable = false
                        end
                    end
                    if (!acceptable)
                        break
                    end
                end
                if (current_error < solution_error)
                    solution_error = current_error
                    problem.camera.position = position
                    problem.camera.matrix = camera_matrix
                end
            end
        end

        for problem in problems
            plot_3dcamera(Plot3dCameraInput(
                problem.camera.euler_rotation,
                problem.camera.position,
            ), :green)
        end

        display("Calculated translation: $(camera_calculated.position)")
        display("Actual translation: $(camera.position)")
        display("Difference between translations: $(translations_difference(camera_calculated.position, camera.position))")

        display("Camera projection matrix: $(camera.matrix ./ camera.matrix[3, 4])")
        display("Calculated projection camera matrix: $(camera_calculated.matrix ./ camera_calculated.matrix[3, 4])")

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

        display(scene.figure)
    end

    function generate_monodromy_solutions()
        scene, problems = create_scene_instances_and_problems(
            random_seed=777,
            number_of_cylinders=4,
            number_of_instances=1,
        )

        display(scene.figure)

        for (i, problem) in enumerate(problems)
            display("--------------- problem $i ---------------")
            for i in 1:4
                begin #asserts
                    lines = problem.lines
                    points_at_infinity = problem.points_at_infinity
                    dualquadrics = problem.dualquadrics
                    R = problem.camera.quaternion_rotation'
                    cameramatrix = problem.camera.matrix
                    top_left_intrinsic = problem.camera.intrinsic[1:3, 1:3]
                    display("$(lines[i, :]' * top_left_intrinsic * R * points_at_infinity[i, :]), (1) line point")
                    display("$(lines[i, :]' * cameramatrix * dualquadrics[i, :, :] * cameramatrix' * lines[i, :]), (2) line quadric")
                end
            end
            display("------------- end problem $i -------------")
        end

        intrinsic_rotation_system = build_intrinsic_rotation_conic_system(problems)
        fₓ, _, _, skew, fᵧ, _, cₓ, cᵧ = vec(problems[1].camera.intrinsic)
        startingsolution = [1, fₓ, fᵧ, skew, cₓ, cᵧ]
        parameters = []
        for problem in problems
            startingsolution = stack_homotopy_parameters(
                startingsolution,
                problem.camera.quaternion_rotation',
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
            max_loops_no_progress = 20,
        )
        display(monodromy_solutions)

        try
            mkdir("tmp")
        catch
        end

        serialize(
            "tmp/intrinsic_rotation_monodromy_solutions.jld",
            MonodromyParametersSolutionsPair(
                parameters,
                solutions(monodromy_solutions)
            )
        )

        translation_system = build_intrinsic_rotation_translation_conic_system(problems[1])
        startingsolution = [1; problems[1].camera.position;]
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
            MonodromyParametersSolutionsPair(
                parameters,
                solutions(monodromy_solutions)
            )
        )
    end

    function create_scene_instances_and_problems(;
        random_seed = 7, 
        number_of_cylinders = 4,
        number_of_instances = 5,
    )
        Random.seed!(random_seed)

        scene = Scene()

        cylinders = []
        for i in 1:number_of_cylinders
            cylinder = Cylinder.CylinderProperties()
            position = normalize(rand(Float64, 3)) * rand_in_range(0.0, 5.0)
            rotation = rand_in_range((0, 180), 3)
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

            standard, dual, singularpoint = Cylinder.standard_and_dual(cylinder.transform, cylinder.radiuses)
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

        instances = []
        intrinsic = build_intrinsic_matrix(IntrinsicParameters(
            focal_length_x = rand_in_range(0.1, 4),
            focal_length_y = rand_in_range(0.1, 4),
            principal_point_x = rand_in_range(-0.5, 0.5),
            principal_point_y = rand_in_range(-0.5, 0.5),
            skew = rand_in_range(0, 1),
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
                conic = Conic.ConicProperties(
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

        plot_3dcylinders(Plot3dCylindersInput(
            [cylinder.transform for cylinder in cylinders],
            [cylinder.radiuses for cylinder in cylinders],
            number_of_cylinders,
        ))

        for (i, instance) in enumerate(instances)
            camera = instance.camera
            conics = instance.conics
            conics_contours = instance.conics_contours
            plot_3dcamera(Plot3dCameraInput(
                camera.euler_rotation,
                camera.position,
            ))
            if i > 1
                add_2d_axis!()
            end
            plot_2dpoints([(conic.singular_point ./ conic.singular_point[3])[1:2] for conic in conics]; axindex = i)
            plot_2dcylinders(conics_contours, alpha=0.5; axindex = i)
        end

        # 3 line minimum to solve the pose
        problems::Vector{CylinderCameraContoursProblem} = []
        numberoflines_tosolvefor = 8
        for instance in instances
            conics_contours = instance.conics_contours

            lines = Matrix{Float64}(undef, numberoflines_tosolvefor, 3)
            points_at_infinity = Matrix{Float64}(undef, numberoflines_tosolvefor, 3)
            dualquadrics = Array{Float64}(undef, numberoflines_tosolvefor, 4, 4)
            possible_picks = 1:(number_of_cylinders*2)
            for store_index in (1:numberoflines_tosolvefor)
                line_index = rand(possible_picks)
                possible_picks = filter(x -> x != line_index, possible_picks)
                i = ceil(Int, line_index / 2)
                j = (line_index - 1) % 2 + 1

                line = conics_contours[i, j, :]
                lines[store_index, :] = normalize(line)
                points_at_infinity[store_index, :] = normalize(cylinders[i].singular_point[1:3])
                dualquadrics[store_index, :, :] = cylinders[i].dual_matrix ./ cylinders[i].dual_matrix[4, 4]
            end

            problem = CylinderCameraContoursProblem(
                instance.camera,
                lines,
                points_at_infinity,
                dualquadrics,
            )
            push!(problems, problem)
        end

        return scene, problems
    end

    export main, generate_monodromy_solutions
end