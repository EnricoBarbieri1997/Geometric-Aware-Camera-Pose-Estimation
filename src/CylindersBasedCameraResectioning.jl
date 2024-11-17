module CylindersBasedCameraResectioning
    include("includes.jl")

    using .Geometry: Line, Cylinder as CylinderType, line_to_homogenous, get_cylinder_contours
    using .Space: transformation, random_transformation, identity_transformation, build_rotation_matrix
    using .Camera: CameraProperties, build_intrinsic_matrix, build_camera_matrix, lookat_rotation
    using .Plotting: initfigure, plot_2dpoints, plot_line_2d, Plot3dCameraInput, plot_3dcamera, Plot3dCylindersInput, plot_cylinders_contours, plot_3dcylinders, plot_2dcylinders
    using .EquationSystems: stack_homotopy_parameters, build_intrinsic_rotation_conic_system, build_intrinsic_rotation_translation_conic_system
    using .Debug
    using .Utils
    using LinearAlgebra: deg2rad, diagm, dot, normalize, svd, pinv
    using HomotopyContinuation, Polynomials, Rotations
    using Random
    using Serialization

    struct MonodromyParametersSolutionsPair
        start_parameters::Vector{Float64}
        solutions::Vector{Vector{ComplexF64}}
    end

    function main()
        Random.seed!(7)

        number_of_cylinders = 4
        cylinders = []
        for i in 1:number_of_cylinders
            cylinder = Cylinder.CylinderProperties()
            position = normalize(rand(Float64, 3)) * rand_in_range(0.0, 5.0)
            rotation = rand_in_range((0, 40), 3)
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

        camera = CameraProperties(
            position = [2.0, 30.0, 5.0],
            euler_rotation = [-83.0, 180.0, 0.0],
        )
        camera.matrix = build_camera_matrix(camera.position, camera.euler_rotation, 2, 1)

        conics = []
        for i in 1:number_of_cylinders
            conic = Conic.ConicProperties(
                pinv(camera.matrix') * cylinders[i].matrix * pinv(camera.matrix),
                camera.matrix * cylinders[i].singular_point,
                camera.matrix * cylinders[i].dual_matrix * camera.matrix',
            )
            push!(conics, conic)
        end

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

        figure = initfigure()
        plot_3dcamera(Plot3dCameraInput(
            camera.euler_rotation,
            camera.position,
        ))
        plot_3dcylinders(Plot3dCylindersInput(
            [cylinder.transform for cylinder in cylinders],
            [cylinder.radiuses for cylinder in cylinders],
            number_of_cylinders,
        ))
        plot_2dpoints([(conic.singular_point ./ conic.singular_point[3])[1:2] for conic in conics])
        # plot_2dcylinders(conics_contours)

        # 3 line minimum to solve the pose
        numberoflines_tosolvefor = 4
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

        parameters_solutions_pair = nothing
        try
            parameters_solutions_pair = deserialize("tmp/intrinsic_rotation_monodromy_solutions.jld")
        catch
            error("generate intrinsic-rotation monodromy first")
            return 1
        end

        intrinsic_rotation_system = build_intrinsic_rotation_conic_system(lines, points_at_infinity)
        parameters = stack_homotopy_parameters(lines, points_at_infinity)

        result = solve(
            intrinsic_rotation_system,
            # parameters_solutions_pair.solutions,
            # start_parameters = parameters_solutions_pair.start_parameters,
            target_parameters = parameters,
        )
        @info result

        camera_calculated = nothing
        solution_error = Inf
        solutions_to_try = [real_solutions(result)[3]]
        for solution in solutions_to_try
            solution = solution ./ solution[5]
            @info solution[5]
            xₛ = solution[1]
            yₛ = solution[2]
            zₛ = solution[3]
            fₛ = solution[4]

            camera_extrinsic_rotation= quat_from_rotmatrix(build_rotation_matrix(xₛ , yₛ, zₛ, true))
            acceptable = true
            current_error = 0
            for (i, contour) in enumerate(eachslice(conics_contours, dims=1))
                for line in eachslice(contour, dims=1)
                    eq = line' * build_intrinsic_matrix(fₛ)[1:3, 1:3] * camera_extrinsic_rotation * cylinders[i].singular_point[1:3]
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
                camera_calculated = CameraProperties(
                    euler_rotation = rad2deg.(Rotations.params(RotXYZ(camera_extrinsic_rotation'))),
                    quaternion_rotation = camera_extrinsic_rotation',
                    focal_length = fₛ,
                )
            end
        end

        begin #asserts
            @assert isnothing(camera_calculated) == false  "(11) Found rotation satisfies constraints"
        end

        display("Error of the best rotation solution: $(solution_error)")
        display("Best solution for rotation: $(camera_calculated.euler_rotation)")
        display("Actual rotation: $(camera.euler_rotation)")

        try
            parameters_solutions_pair = deserialize("tmp/translation_monodromy_solutions.jld")
        catch
            error("generate translation monodromy first")
            return 1
        end

        translation_system = build_intrinsic_rotation_translation_conic_system(
            build_intrinsic_matrix(camera_calculated.focal_length),
            camera_calculated.quaternion_rotation,
            lines[1:3, :],
            dualquadrics[1:3, :, :]
        )

        parameters = stack_homotopy_parameters(lines[1:3, :], dualquadrics[1:3, :, :])

        result = solve(
            translation_system,
            # parameters_solutions_pair.solutions,
            # start_parameters = parameters_solutions_pair.start_parameters,
            target_parameters = parameters,
        )
        @info result

        solution_error = Inf
        for solution in real_solutions(result)
            solution = solution ./ solution[4]
            txₛ = solution[1]
            tyₛ = solution[2]
            tzₛ = solution[3]
            position = [txₛ, tyₛ, tzₛ]
            acceptable = true
            current_error = 0
            Pₛ = build_camera_matrix(position, camera_calculated.quaternion_rotation, camera_calculated.focal_length, 1)
            for (i, conic_contour) in enumerate(eachslice(conics_contours, dims=1))
                for line in eachslice(conic_contour, dims=1)
                    eq = line' * Pₛ * cylinders[i].dual_matrix * Pₛ' * line
                    current_error += abs(eq)
                    # if (!(eq ≃ 0))
                    #     acceptable = false
                    # end
                    if (!acceptable)
                        break
                    end
                end
            end
            if (current_error < solution_error)
                solution_error = current_error
                camera_calculated.position = position
            end
        end

        begin #asserts
            @assert solution_error < Inf  "(12) Found translation satisfies constraints"
        end

        camera_calculated.matrix = build_camera_matrix(camera_calculated.position, camera_calculated.quaternion_rotation, camera_calculated.focal_length, 1)

        @info "Calculated translation: $(camera_calculated.position)"
        @info "Actual translation: $(camera.position)"

        display("Camera projection matrix: $(camera.matrix ./ camera.matrix[3, 4])")
        display("Calculated projection camera matrix: $(camera_calculated.matrix ./ camera_calculated.matrix[3, 4])")

        reconstructued_contours = Array{Float64}(undef, number_of_cylinders, 2, 3)
        for i in 1:number_of_cylinders
            lines = get_cylinder_contours(
                cylinders[i].geometry,
                collect(camera_calculated.position),
                camera_calculated.matrix
            )
            for (j, line) in enumerate(lines)
                line_homogenous = line_to_homogenous(line)
                reconstructued_contours[i, j, :] = line_homogenous
            end
        end

        plot_2dcylinders(reconstructued_contours, linestyle=:dash)

        figure
    end

    function generate_monodromy_solutions()
        Random.seed!(777)

        camera_translationdirection = normalize(rand(Float64, 3))
        camera_translation = camera_translationdirection * rand_in_range(15.0, 20.0)
        # camera_rotation = lookat_rotation(camera_translationdirection, [0, 0, 0], [0, 0, 1])
        rx, ry, rz = (-90, -110, 0)
        camera_rotation = RotXYZ(deg2rad(rx), deg2rad(ry), deg2rad(rz))

        # quaternion_camera_rotation = quat_from_rotmatrix(camera_rotation)
        quaternion_camera_rotation = QuatRotation(camera_rotation)
        rotation_params = Rotations.params(quaternion_camera_rotation')
        rotation_params = rotation_params ./ rotation_params[1]
        w, x, y, z = rotation_params
        f = 2

        cameramatrix = build_camera_matrix(
            camera_translation,
            quaternion_camera_rotation, f, 1
        )

        random_cylindertranslation_range = 5
        cylinders = []
        cylinders_geometry = []
        transforms::Vector{Matrix{Float64}} = []
        radiuses::Vector{Vector{Float64}} = []
        for i in 1:2
            position_direction = normalize(rand(Float64, 3))
            position = position_direction * rand_in_range(0, random_cylindertranslation_range)
            rotation = rand_in_range((0, 90), 3)

            transform = transformation(position, rotation)
            push!(transforms, transform)
            radius = rand_in_range(1, 3)
            radius = [radius, radius]
            push!(radiuses, radius)

            axis = transform * [0; 0; 1; 0]
            axis = axis ./ axis[3]
            axis = axis[1:3]
            push!(cylinders_geometry, CylinderType(
                position,
                radius[1],
                axis,
            ))

            push!(cylinders, Cylinder.standard_and_dual(transform, radius))
        end

        contour_lines = Array{Float64}(undef, 2, 2, 3)
        lines = Matrix{Float64}(undef, 4, 3)
        points_at_infinity = Matrix{Float64}(undef, 4, 3)
        dualquadrics = Array{Float64}(undef, 4, 4, 4)
        for i in 1:2
            cylinder_lines = get_cylinder_contours(
                cylinders_geometry[i],
                camera_translation,
                cameramatrix
            )
            for (line_index, line) in enumerate(cylinder_lines)
                line_homogenous = line_to_homogenous(line)
                index = (i - 1) * 2 + line_index
                contour_lines[i, line_index, :] = line_homogenous
                lines[index, :] = normalize(line_homogenous)
                points_at_infinity[index, :] = normalize(cylinders[i][3][1:3])
                dualquadrics[index, :, :] = cylinders[i][2] ./ cylinders[i][2][4, 4]
            end
        end

        figure = initfigure()
        # rx, ry, rz = rad2deg.(Rotations.params(RotXYZ(quaternion_camera_rotation)))
        plot_3dcamera(Plot3dCameraInput(
            [rx, ry, rz],
            camera_translation,
        ))
        plot_3dcylinders(Plot3dCylindersInput(
            transforms,
            radiuses,
            2,
        ))
        plot_2dcylinders(contour_lines)

        plot_2dpoints([
            cameramatrix * [1, 0, 0, 1],
            cameramatrix * [0, 1, 0, 1],
            cameramatrix * [0, 0, 1, 1],
        ])

        display(figure)

        intrinsic_rotation_system = build_intrinsic_rotation_conic_system(lines, points_at_infinity)
        parameters = stack_homotopy_parameters(lines, points_at_infinity)
        startingsolution = [x, y, z, f]
        monodromy_solutions = monodromy_solve(intrinsic_rotation_system, startingsolution, parameters, max_loops_no_progress=5)
        @info monodromy_solutions

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

        intrinsic_matrix = build_intrinsic_matrix(f)

        translation_system = build_intrinsic_rotation_translation_conic_system(
            intrinsic_matrix,
            quaternion_camera_rotation,
            lines[1:3, :],
            dualquadrics[1:3, :, :]
        )
        parameters = stack_homotopy_parameters(lines[1:3, :], dualquadrics[1:3, :, :])
        monodromy_solutions = monodromy_solve(translation_system, camera_translation, parameters, max_loops_no_progress=5)
        @info monodromy_solutions

        serialize(
            "tmp/translation_monodromy_solutions.jld",
            MonodromyParametersSolutionsPair(
                parameters,
                solutions(monodromy_solutions)
            )
        )
    end

    export main, generate_monodromy_solutions
end
