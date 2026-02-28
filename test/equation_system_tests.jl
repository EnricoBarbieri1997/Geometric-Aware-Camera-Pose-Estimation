using CylindersBasedCameraResectioning.Camera: CameraProperties, CameraViewPair, random_camera_lookingat_center
using CylindersBasedCameraResectioning.Geometry: get_view
using CylindersBasedCameraResectioning.Cylinder: CylinderProperties, CalibrationRigs, points_at_infinity_dualquadrics
using CylindersBasedCameraResectioning.IO: read_axis_rig_lines, read_camera
using CylindersBasedCameraResectioning.EquationSystems: build_intrinsic_rotation_conic_system, stack_homotopy_parameters
using CylindersBasedCameraResectioning.EquationSystems.Minimization: build_intrinsic_rotation_conic_system as build_intrinsic_rotation_conic_system_minimization
using CylindersBasedCameraResectioning.EquationSystems.Problems: CylinderCameraContoursProblem, CylinderCameraContoursProblemValidationData
using CylindersBasedCameraResectioning.EquationSystems.Problems.IntrinsicParameters: Configurations as IntrinsicParametersConfigurations
using CylindersBasedCameraResectioning.Scene: intrinsic_rotation_system_setup
using CylindersBasedCameraResectioning.Utils: rad2deg, eulerangles_from_rotationmatrix, lines_clp_to_stack

using HomotopyContinuation
using LinearAlgebra: norm
using Rotations



@testset "solution_solve_system" begin
    intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ_cₓ_cᵧ
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)

    intrinsics = [
        2666.6666666666665 0.0 960.0;
        0.0 2666.6666666666665 540.0;
        0.0 0.0 1.0
    ]
    
    cameras = []
    for i in 1:2
        position, rotation = random_camera_lookingat_center()
        camera = CameraProperties()
        camera.position = position
        camera.quaternion_rotation = rotation
        camera.intrinsic = intrinsics
        push!(cameras, camera)
    end
    camera1 = cameras[1]
    camera2 = cameras[2]

    view_1 = get_view(cylinders, camera1)
    view_2 = get_view(cylinders, camera2)

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
    )
    for (i, camera_view_pair) in enumerate(camera_view_pairs)
        lines = lines_clp_to_stack(camera_view_pair.view)
        problem = CylinderCameraContoursProblem(
            camera_view_pair.camera,
            lines,
            lines,
            points_at_infinity,
            dualquadrics,
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
            lines',
        )
    end
    parameters = convert(Vector{Float64}, parameters)

    rot1 = Rotations.params(QuatRotation(cameras[1].rotation_matrix))
    rot1 = rot1 / rot1[1]
    rot1 = rot1[2:4]
    rot2 = Rotations.params(QuatRotation(cameras[2].rotation_matrix))
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
@testset "solution_solve_minimization_system" begin
    intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ_cₓ_cᵧ
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)

    intrinsics = [
        2666.6666666666665 0.0 960.0;
        0.0 2666.6666666666665 540.0;
        0.0 0.0 1.0
    ]
    
    cameras = []
    for i in 1:2
        position, rotation = random_camera_lookingat_center()
        camera = CameraProperties()
        camera.position = position
        camera.quaternion_rotation = rotation
        camera.intrinsic = intrinsics
        push!(cameras, camera)
    end
    camera1 = cameras[1]
    camera2 = cameras[2]

    view_1 = get_view(cylinders, camera1)
    view_2 = get_view(cylinders, camera2)

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
    )
    for (i, camera_view_pair) in enumerate(camera_view_pairs)
        lines = lines_clp_to_stack(camera_view_pair.view)
        problem = CylinderCameraContoursProblem(
            camera_view_pair.camera,
            lines,
            lines,
            points_at_infinity,
            dualquadrics,
            validation_data,
            UInt8(intrinsic_configuration)
        )
        push!(problems, problem)
    end

    rotation_intrinsic_system = build_intrinsic_rotation_conic_system_minimization(
        problems;
    )
    parameters = []
    for problem in problems
        lines = problem.lines[1:5, :]
        parameters = stack_homotopy_parameters(
            parameters,
            lines',
        )
    end
    parameters = convert(Vector{Float64}, parameters)

    rot1 = Rotations.params(QuatRotation(cameras[1].rotation_matrix))
    rot1 = rot1 / rot1[1]
    rot1 = rot1[2:4]
    rot2 = Rotations.params(QuatRotation(cameras[2].rotation_matrix))
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

function solution_solve_system_with_auto_params_base(; minimization=false)
    intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ_cₓ_cᵧ
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)

    intrinsics = [
        2666.6666666666665 0.0 960.0;
        0.0 2666.6666666666665 540.0;
        0.0 0.0 1.0
    ]
    
    cameras = []
    for i in 1:2
        position, rotation = random_camera_lookingat_center()
        camera = CameraProperties()
        camera.position = position
        camera.quaternion_rotation = rotation
        camera.intrinsic = intrinsics
        push!(cameras, camera)
    end
    camera1 = cameras[1]
    camera2 = cameras[2]

    view_1 = get_view(cylinders, camera1)
    view_2 = get_view(cylinders, camera2)

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
    )
    for (i, camera_view_pair) in enumerate(camera_view_pairs)
        lines = lines_clp_to_stack(camera_view_pair.view)
        problem = CylinderCameraContoursProblem(
            camera_view_pair.camera,
            lines,
            lines,
            points_at_infinity,
            dualquadrics,
            validation_data,
            UInt8(intrinsic_configuration)
        )
        push!(problems, problem)
    end

    rotation_intrinsic_system, parameters = intrinsic_rotation_system_setup(
        problems;
        minimization=false,
        intrinsic_configuration
    )

    rot1 = Rotations.params(QuatRotation(cameras[1].rotation_matrix))
    rot1 = rot1 / rot1[1]
    rot1 = rot1[2:4]
    rot2 = Rotations.params(QuatRotation(cameras[2].rotation_matrix))
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
@testset "solution_solve_system_with_auto_params" begin
    solution_solve_system_with_auto_params_base()
end
@testset "solution_solve_minimization_system_with_auto_params" begin
    solution_solve_system_with_auto_params_base(; minimization=true)
end

function solution_solve_system_with_auto_params_with_combination(; minimization=false)
    intrinsic_configuration = IntrinsicParametersConfigurations.fₓ_fᵧ_cₓ_cᵧ
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)

    intrinsics = [
        2666.6666666666665 0.0 960.0;
        0.0 2666.6666666666665 540.0;
        0.0 0.0 1.0
    ]
    
    cameras = []
    for i in 1:2
        position, rotation = random_camera_lookingat_center()
        camera = CameraProperties()
        camera.position = position
        camera.quaternion_rotation = rotation
        camera.intrinsic = intrinsics
        push!(cameras, camera)
    end
    camera1 = cameras[1]
    camera2 = cameras[2]

    view_1 = get_view(cylinders, camera1)
    view_2 = get_view(cylinders, camera2)

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
    )
    for (i, camera_view_pair) in enumerate(camera_view_pairs)
        lines = lines_clp_to_stack(camera_view_pair.view)
        problem = CylinderCameraContoursProblem(
            camera_view_pair.camera,
            lines,
            lines,
            points_at_infinity,
            dualquadrics,
            validation_data,
            UInt8(intrinsic_configuration)
        )
        push!(problems, problem)
    end

    combination = [
        3, 7, 4, 12, 5,
        8, 6, 9, 1, 10,
    ]

    rotation_intrinsic_system, parameters = intrinsic_rotation_system_setup(
        problems;
        minimization=false,
        intrinsic_configuration,
        equation_combinations=combination,
    )

    rot1 = Rotations.params(QuatRotation(cameras[1].rotation_matrix))
    rot1 = rot1 / rot1[1]
    rot1 = rot1[2:4]
    rot2 = Rotations.params(QuatRotation(cameras[2].rotation_matrix))
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
@testset "solution_solve_system_with_auto_params_with_combination" begin
    solution_solve_system_with_auto_params_with_combination()
end
@testset "solution_solve_minimization_system_with_auto_params_with_combination" begin
    solution_solve_system_with_auto_params_with_combination(; minimization=true)
end