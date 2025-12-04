module SolutionGenerator
    export generate_configurations, generate_parameter_solution_pair, generate_all
    using ..Geometry: get_cylinder_contours
    using ..Space: transformation
    using ..Cylinder: CylinderProperties, standard_and_dual as standard_and_dual_cylinder
    using ..Camera: CameraProperties, lookat_rotation
    using ..Scene: SceneData, InstanceConfiguration, intrinsic_rotation_system_setup
    using ..EquationSystems.Problems: CylinderCameraContoursProblem, CylinderCameraContoursProblemValidationData
    using ..EquationSystems.Problems.IntrinsicParameters: Configurations as IntrinsicParametersConfigurations
    using ..Utils: eulerangles_from_rotationmatrix
    using LinearAlgebra, Serialization 
    using HomotopyContinuation, Rotations

    ASSERTS_ENABLED = false

    function generate_configurations()
        quadrants_directions::Array{Vector{Float64}} = [
            normalize([1.0, 1.0, 1.0]),
            normalize([-1.0, 1.0, 1.0]),
            normalize([-1.0, 1.0, -1.0]),
            normalize([-1.0, -1.0, 1.0]),
            normalize([-1.0, -1.0, -1.0]),
            normalize([1.0, -1.0, 1.0]),
            normalize([1.0, -1.0, -1.0]),
            normalize([1.0, 1.0, -1.0]),
        ]

        cameras = []
        for dir in quadrants_directions
            camera = CameraProperties()
            camera.position = dir * 10.0
            rotation_matrix = lookat_rotation(camera.position, [0.0, 0.0, 0.0])
            camera.quaternion_rotation = QuatRotation(rotation_matrix)
			camera.euler_rotation = rad2deg.(eulerangles_from_rotationmatrix(rotation_matrix))
            camera.intrinsic = [
                2666.6666 0 960;
                0 1500 540;
                0 0 1
            ]
            push!(cameras, camera)
        end

        cylinders = canonical_cylinder_rig()

        views::Array{Array{Float64, 3}} = []
        for camera in cameras
            conics_contours = Array{Float64}(undef, 3, 2, 3)
            for (i, cylinder) in enumerate(cylinders)
                lines = get_cylinder_contours(
                    cylinder,
                    camera
                )
                for (j, line) in enumerate(lines)
                    conics_contours[i, j, :] = line

                    if (ASSERTS_ENABLED)
                        @assert line' * conics[i].dual_matrix * line ≃ 0 "(3) Line of projected singular plane $(1) belongs to the dual conic $(1)"
                        @assert line' * camera.matrix * cylinders[i].singular_point ≃ 0 "(8) Line $(j) of conic $(i) passes through the projected singular point"
                        err = (line' * camera.matrix * cylinders[i].dual_matrix * camera.matrix' * line)
                        @assert err ≃ 0 "(9) Line $(j) of conic $(i) is tangent to the projected cylinder. $(err)"
                    end
                end
            end
            push!(views, conics_contours)
        end

        for (i, camera) in enumerate(cameras)
            camera_view_pair = CameraViewPair(i, camera, views[i])
            serialize("./tmp/start_solutions/camera_view_pairs/$(i).jls", camera_view_pair)
        end
    end

    function generate_parameter_solution_pair(index::String)
        intrinsic_configuration = IntrinsicParametersConfigurations.focal_length_x | IntrinsicParametersConfigurations.focal_length_y | IntrinsicParametersConfigurations.principal_point_x | IntrinsicParametersConfigurations.principal_point_y

        cylinders = canonical_cylinder_rig()
        camera_index_1, camera_index_2 = parse.(Int, split(index, "_"))
        camera_view_pair_1 = deserialize("./tmp/start_solutions/camera_view_pairs/$(camera_index_1).jls")
        camera_view_pair_2 = deserialize("./tmp/start_solutions/camera_view_pairs/$(camera_index_2).jls")

        scene = SceneData()
        scene.cylinders = cylinders
        
        scene.instances = map(camera_view_pair -> begin
            instance = InstanceConfiguration()
            instance.camera = camera_view_pair.camera
            instance.conics_contours = camera_view_pair.view
            return instance
        end, [camera_view_pair_1, camera_view_pair_2])

        points_at_infinity::Matrix{Float64} = zeros(Float64, 6, 3)
        for (i, cylinder) in enumerate(cylinders)
            doubled_index = (i-1)*2+1
            points_at_infinity[doubled_index, :] = normalize(cylinder.singular_point[1:3])
            points_at_infinity[doubled_index+1, :] = normalize(cylinder.singular_point[1:3])
        end
        dualquadrics::Array{Float64, 3} = zeros(Float64, 6, 4, 4)
        for (i, cylinder) in enumerate(cylinders)
            doubled_index = (i-1)*2+1
            dualquadrics[doubled_index, :, :] = cylinder.dual_matrix ./ cylinder.dual_matrix[4, 4]
            dualquadrics[doubled_index+1, :, :] = cylinder.dual_matrix ./ cylinder.dual_matrix[4, 4]
        end

        problems::Vector{CylinderCameraContoursProblem} = []
        validation_data = CylinderCameraContoursProblemValidationData(
            Matrix{Float64}(undef, 0, 3),
            Matrix{Float64}(undef, 0, 3),
            Array{Float64}(undef, 0, 4, 4),
            Vector{Float64}(undef, 6)
        )
        for (i, camera_view_pair) in enumerate([camera_view_pair_1, camera_view_pair_2])
            lines = reshape(camera_view_pair.view, 6, 3)
            if (i == 2)
                lines = lines[1:4, :]
            end
            problem = CylinderCameraContoursProblem(
                camera_view_pair.camera,
                lines,
                lines,
                points_at_infinity,
                dualquadrics,
                collect(1:size(lines)[1]),
                validation_data,
                intrinsic_configuration
            )
            push!(problems, problem)
        end

        system, parameters = intrinsic_rotation_system_setup(problems; intrinsic_configuration)
        result = solve(
            system;
            target_parameters=parameters,
            start_system = :total_degree,
        )

        sol = solutions(result)
        if length(sol) < 128
            display("Warning: The number of solutions found is less than expected: $(length(sol)) < 128 for index $(index)")
        end
        parameters_solutions_pair = ParameterSolutionsPair(
            index,
            parameters,
            sol
        )

        serialize("./tmp/start_solutions/parameters_solution_pairs/$(index).jls", parameters_solutions_pair)
    end

    function generate_all()
        generate_configurations()
        for i in 1:7
            for j in i+1:8
                if (i == 1 && j == 2) # Only during testing because we already have this one
                    continue
                end
                index = "$(i)_$(j)"
                generate_parameter_solution_pair(index)
            end
        end
    end

    struct CameraViewPair
        index::Int
        camera::CameraProperties
        view::Array{Float64, 3}
    end

    struct ParameterSolutionsPair
        index::String
        parameters::Vector{Float64}
        solutions::Vector{Vector{Number}}
    end

    function rig_cylinder(euler_rotation)
        cylinder = CylinderProperties()
        position = [0.0, 0.0, 0.0]
        cylinder.euler_rotation = euler_rotation
        cylinder.transform = transformation(position, cylinder.euler_rotation)
        radius = [1.0, 1.0]
        cylinder.radiuses = [radius[1], radius[1]]

        axis = cylinder.transform * [0; 0; 1; 0]
        axis = axis ./ axis[3]
        axis = axis[1:3]

        standard, dual, singularpoint = standard_and_dual_cylinder(cylinder.transform, cylinder.radiuses)
        cylinder.matrix = standard
        cylinder.dual_matrix = dual
        cylinder.singular_point = singularpoint

        if (ASSERTS_ENABLED)
            @assert cylinder.matrix * cylinder.dual_matrix ≃ diagm([1, 1, 0, 1]) "(-1) The dual quadric is correct"

            @assert cylinder.singular_point' * cylinder.matrix * cylinder.singular_point ≃ 0 "(1) Singular point $(1) belongs to the cylinder $(1)"
            dual_singular_plane = cylinder.transform * reshape([1, 0, 0, -cylinder.radiuses[1]], :, 1)
            @assert (dual_singular_plane' * cylinder.dual_matrix * dual_singular_plane) ≃ 0 "(2) Perpendicular plane $(1) belongs to the dual cylinder $(1)"

            @assert (cylinder.matrix * cylinder.singular_point) ≃ [0, 0, 0, 0] "(6) Singular point is right null space of cylinder matrix $(i)"

            @assert ((dual_singular_plane' * cylinder.singular_point) ≃ 0 && (dual_singular_plane' * cylinder.dual_matrix * dual_singular_plane) ≃ 0) "(7) Singular plane / point and dual quadric constraints $(i)"
            @assert cylinder.singular_point[4] ≃ 0 "(10) Singular point is at infinity $(i)"
        end

        return cylinder
    end

    function canonical_cylinder_rig()
        cylinders = []
        # Cylinder X
        push!(cylinders, rig_cylinder([0.0, 90.0, 0.0]))
        # Cylinder Y
        push!(cylinders, rig_cylinder([90.0, 0.0, 0.0]))
        # Cylinder Z
        push!(cylinders, rig_cylinder([0.0, 0.0, 0.0]))
        return cylinders
    end

end