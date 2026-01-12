module SolutionGenerator
export generate_configurations, generate_parameter_solution_pair, generate_all
using ..Geometry: get_view, homogeneous_anglebetween
using ..Cylinder: CylinderProperties, CalibrationRigs, points_at_infinity_dualquadrics
using ..Camera: CameraProperties, CameraViewPair, lookat_rotation
using ..Scene: SceneData, InstanceConfiguration, intrinsic_rotation_system_setup, create_scene_instances_and_problems, plot_reconstructed_scene, plot_scene
using ..EquationSystems.Problems: CylinderCameraContoursProblem, CylinderCameraContoursProblemValidationData
using ..EquationSystems.Problems.IntrinsicParameters: Configurations as IntrinsicParametersConfigurations
using ..Utils: eulerangles_from_rotationmatrix, rand_in_range
using ..Plotting: initfigure, plot_2dcylinders
using LinearAlgebra, Serialization
using HomotopyContinuation, Rotations
using Random

CAMERAS_POOL_SIZE = 8

function solve_by_similarity()
    Random.seed!(785687)
    intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ_cₓ_cᵧ
    cylinders::Vector{CylinderProperties} = CalibrationRigs.axis_rig()

    intrinsics = [
        rand_in_range(2500.0, 2700.0) 0.0 rand_in_range(950.0, 970.0);   # fₓ, skew, cₓ
        0.0 rand_in_range(1400.0, 1600.0) rand_in_range(530.0, 550.0);   # 0, fᵧ, cᵧ
        0.0 0.0 1.0                                                      # bottom row
    ]

    cameras::Vector{CameraProperties} = []
    for i in 1:2
        dir = normalize(Random.randn(3))
        camera = CameraProperties()
        camera.position = dir * rand_in_range(10.0, 15.0)
        # camera.position = [10.0, 10.0, 0.0]
        # if i == 2
        #     camera.position = camera.position * -1
        # end
        rotation_matrix = lookat_rotation(camera.position, [0.0, 0.0, 0.0])
        camera.quaternion_rotation = QuatRotation(rotation_matrix)
        camera.euler_rotation = rad2deg.(eulerangles_from_rotationmatrix(rotation_matrix))
        camera.intrinsic = intrinsics
        push!(cameras, camera)
    end

    views::Array{Array{Float64,3}} = []
    for camera in cameras
        conics_contours = get_view(cylinders, camera)
        push!(views, conics_contours)
    end

    camera_view_pairs = [
        CameraViewPair(1, cameras[1], views[1]),
        CameraViewPair(2, cameras[2], views[2]),
    ]

    scene = SceneData()
    scene.figure = initfigure()
    scene.cylinders = cylinders

    scene.instances = map(camera_view_pair -> begin
            instance = InstanceConfiguration()
            instance.camera = camera_view_pair.camera
            instance.conics_contours = camera_view_pair.view
            return instance
        end, camera_view_pairs)

    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)

    original_problems::Vector{CylinderCameraContoursProblem} = []
    validation_data = CylinderCameraContoursProblemValidationData(
        Matrix{Float64}(undef, 0, 3),
        Matrix{Float64}(undef, 0, 3),
        Array{Float64}(undef, 0, 4, 4),
        Vector{Float64}(undef, 6)
    )
    for (i, camera_view_pair) in enumerate(camera_view_pairs)
        lines = reshape(camera_view_pair.view, 6, 3)
        # if (i == 1)
        #     lines = lines[1:4, :]
        # end
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
            UInt8(intrinsic_configuration)
        )
        push!(original_problems, problem)
    end

    plot_scene(scene, original_problems)

    rotation_intrinsic_system, parameters = intrinsic_rotation_system_setup(
        original_problems;
        minimization=false,
        intrinsic_configuration
    )

    rot1 = Rotations.params(scene.instances[1].camera.quaternion_rotation)
    rot1 = rot1 / rot1[1]
    rot1 = rot1[2:4]
    rot2 = Rotations.params(scene.instances[2].camera.quaternion_rotation)
    rot2 = rot2 / rot2[1]
    rot2 = rot2[2:4]

    intrinsics = scene.instances[1].camera.intrinsic
    factor = 1.0 # / 3000.0
    solution = [
        intrinsics[1, 1] * factor,
        intrinsics[2, 2] * factor,
        intrinsics[1, 3] * factor,
        intrinsics[2, 3] * factor,
        rot1...,
        rot2...,
    ]
    # display(pair[3:4])
    display("Known solution")
    display(solution;)

    eq = expressions(rotation_intrinsic_system)[1]
    display(eq)
    display(original_problems[1].lines[1, :])

    # display(evaluate(rotation_intrinsic_system, solution, parameters))

    # display(scene.figure)

    return

    pairs = camera_pairs_by_similarity(scene.instances[1].conics_contours, scene.instances[2].conics_contours)

    for pair in pairs
        reference_start = deserialize("./tmp/start_solutions/parameters_solution_pairs/$(pair[1]).jls")
        problems = original_problems[pair[3:4]]

        ids = parse.(Int, split(pair[1], "_"))
        cv_1 = deserialize("./tmp/start_solutions/camera_view_pairs/$(ids[1]).jls")
        cv_2 = deserialize("./tmp/start_solutions/camera_view_pairs/$(ids[2]).jls")

        # plot_2dcylinders(cv_1.view; axindex = 1, linestyle = :dash,)
        # plot_2dcylinders(cv_2.view; axindex = 2, linestyle = :dash,)

        rotation_intrinsic_system, parameters = intrinsic_rotation_system_setup(
            problems;
            minimization=false,
            intrinsic_configuration
        )

        # display("Total degree solution")
        # tdsol = solve(
        #     rotation_intrinsic_system;
        #     start_system = :total_degree,
        #     target_parameters = parameters,
        # )
        # display(tdsol)

        result = solve(
            rotation_intrinsic_system,
            reference_start.solutions;
            start_parameters=reference_start.parameters,
            target_parameters=parameters,
        )

        @info result
    end

    display(scene.figure)
end

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

    cylinders = CalibrationRigs.axis_rig()

    views::Array{Array{Float64,3}} = []
    for camera in cameras
        conics_contours = get_view(cylinders, camera)
        push!(views, conics_contours)
    end

    for (i, camera) in enumerate(cameras)
        camera_view_pair = CameraViewPair(i, camera, views[i])
        serialize("./tmp/start_solutions/camera_view_pairs/$(i).jls", camera_view_pair)
    end
end

function generate_parameter_solution_pair(index::String)
    intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ_cₓ_cᵧ

    cylinders = CalibrationRigs.axis_rig()
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

    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)

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
            UInt8(intrinsic_configuration)
        )
        push!(problems, problem)
    end

    system, parameters = intrinsic_rotation_system_setup(problems; intrinsic_configuration)
    result = solve(
        system;
        target_parameters=parameters,
        start_system=:total_degree,
    )

    sol = solutions(result)
    # additional_solutions = solutions(monodromy_solve(system, sol, parameters))
    # sol = vcat(sol, additional_solutions)
    if length(sol) < 128
        display("Warning: The number of solutions found is less than expected: $(length(sol)) < 128 for index $(index)")
    end
    parameters_solutions_pair = ParameterSolutionsPair(
        index,
        parameters,
        sol
    )

    display(system)
    display(evaluate(system, sol[1], parameters))

    serialize("./tmp/start_solutions/parameters_solution_pairs/$(index).jls", parameters_solutions_pair)
end

function generate_all()
    generate_configurations()
    for i in 1:(CAMERAS_POOL_SIZE-1)
        for j in i+1:CAMERAS_POOL_SIZE
            if (i == 1 && j == 2) # Only during testing because we already have this one
                continue
            end
            index = "$(i)_$(j)"
            generate_parameter_solution_pair(index)
        end
    end
end

function camera_indices_by_similarity(view)
    camera_scores = fill(typemax(Float64), CAMERAS_POOL_SIZE)
    for camera_index in 1:CAMERAS_POOL_SIZE
        try
            camera_view_pair = deserialize("./tmp/start_solutions/camera_view_pairs/$(camera_index).jls")
            view_stored = camera_view_pair.view
            score = 0.0
            for i in 1:3
                for j in 1:2
                    line_stored = view_stored[i, j, :]
                    line_target = view[i, j, :]

                    score += abs(homogeneous_anglebetween(line_stored, line_target))
                end
            end
            camera_scores[camera_index] = score
        catch e
            display("Error in camera $(camera_index)")
            Base.showerror(stdout, e)
            Base.show_backtrace(stdout, catch_backtrace())
        end
    end
    sorted_indices = sortperm(camera_scores, rev=true)
    merged = [
        [sorted_index, camera_scores[i]] for (i, sorted_index) in enumerate(sorted_indices)
    ]
    return merged
end
function camera_pairs_by_similarity(view_1, view_2)
    indices_1 = camera_indices_by_similarity(view_1)
    indices_2 = camera_indices_by_similarity(view_2)

    done_indices = []
    pairs = []

    for index_1 in indices_1
        for index_2 in indices_2
            if index_1[1] == index_2[1]
                continue
            end
            pair = [Int64(index_1[1]), Int64(index_2[1])]
            indices_to_pick = sort(pair)
            view_pair = pair[1] == indices_to_pick[1] ? [1, 2] : [2, 1]
            score = index_1[2] + index_2[2]
            index = join(indices_to_pick, "_")

            if !(index in done_indices)
                push!(pairs, [index, score, view_pair[1], view_pair[2]])
                push!(done_indices, index)
            end
        end
    end

    return sort(pairs, by=x -> x[2])
end

struct ParameterSolutionsPair
    index::String
    parameters::Vector{Float64}
    solutions::Vector{Vector{Number}}
end
end