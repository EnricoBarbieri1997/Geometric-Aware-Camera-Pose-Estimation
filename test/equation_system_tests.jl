using CylindersBasedCameraResectioning.Camera: CameraViewPair
using CylindersBasedCameraResectioning.Cylinder: CylinderProperties, CalibrationRigs, points_at_infinity_dualquadrics
using CylindersBasedCameraResectioning.IO: read_axis_rig_lines, read_camera
using CylindersBasedCameraResectioning.EquationSystems: build_intrinsic_rotation_conic_system, stack_homotopy_parameters
using CylindersBasedCameraResectioning.EquationSystems.Problems: CylinderCameraContoursProblem, CylinderCameraContoursProblemValidationData
using CylindersBasedCameraResectioning.EquationSystems.Problems.IntrinsicParameters: Configurations as IntrinsicParametersConfigurations

using HomotopyContinuation
using LinearAlgebra: norm
using Rotations

@testset "do_solution_solve_system" begin
    intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ_cₓ_cᵧ
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)
    camera1 = read_camera("../assets/test_scenes/axis_rig/scene.json"; object_path="0.camera")
    camera2 = read_camera("../assets/test_scenes/axis_rig/scene.json"; object_path="1.camera")
    cameras = [camera1, camera2]
    lines_view_1 = read_axis_rig_lines("../assets/test_scenes/axis_rig/scene.json"; object_path="0.contours")
    lines_view_2 = read_axis_rig_lines("../assets/test_scenes/axis_rig/scene.json"; object_path="1.contours")

    view_1 = Array{Float64,3}(undef, 3, 2, 3)
    view_2 = Array{Float64,3}(undef, 3, 2, 3)

    for i in 1:3
        for j in 1:2
            view_1[i, j, :] = lines_view_1[i][j][1:3]
            view_2[i, j, :] = lines_view_2[i][j][1:3]
        end
    end


    camera_view_pairs = [
        CameraViewPair(
            1,
            camera1,
            view_1
        ),
        CameraViewPair(
            2,
            camera2,
            view_2
        )
    ]
    problems::Vector{CylinderCameraContoursProblem} = []
    validation_data = CylinderCameraContoursProblemValidationData(
        Matrix{Float64}(undef, 0, 3),
        Matrix{Float64}(undef, 0, 3),
        Array{Float64}(undef, 0, 4, 4),
        Vector{Float64}(undef, 6)
    )
    for (i, camera_view_pair) in enumerate(camera_view_pairs)
        lines = reshape(camera_view_pair.view, 6, 3)
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

    rotation_intrinsic_system = build_intrinsic_rotation_conic_system(
        problems;
    )
    parameters = []
    for problem in problems
        lines = problem.lines[1:5, :]
        parameters = stack_homotopy_parameters(
            parameters,
            lines,
        )
    end
    parameters = convert(Vector{Float64}, parameters)

    rot1 = Rotations.params(cameras[1].rotation_matrix)
    rot1 = rot1 / rot1[1]
    rot1 = rot1[2:4]
    rot2 = Rotations.params(cameras[2].rotation_matrix)
    rot2 = rot2 / rot2[1]
    rot2 = rot2[2:4]

    intrinsics = cameras[1].intrinsic
    factor = 1.0 / 3000.0
    solution = [
        intrinsics[1, 1] * factor,
        intrinsics[2, 2] * factor,
        intrinsics[1, 3] * factor,
        intrinsics[2, 3] * factor,
        rot1...,
        rot2...,
    ]

    equation_results = evaluate(rotation_intrinsic_system, solution, parameters)
    @test isapprox(norm(equation_results), 0.0; atol=1e-6)
end