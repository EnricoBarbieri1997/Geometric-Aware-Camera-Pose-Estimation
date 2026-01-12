using CylindersBasedCameraResectioning.Space: build_rotation_matrix
using CylindersBasedCameraResectioning.Cylinder: CalibrationRigs, points_at_infinity_dualquadrics, get_view
using CylindersBasedCameraResectioning.Camera: CameraProperties, lookat_rotation, random_camera_lookingat_center
using CylindersBasedCameraResectioning.IO: read_axis_rig_lines, read_camera

using LinearAlgebra: I
using Rotations

intrinsic = [
    2666.6666666666665 0.0 960.0;
    0.0 2666.6666666666665 540.0;
    0.0 0.0 1.0
]

@testset "unrotated_calibrated" begin
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)
    position = [10.0, 10.0, -10.0]
    rotation_matrix = Matrix{Float64}(I, 3, 3)
    quaternion_camera_rotation = QuatRotation(rotation_matrix)
    euler_rotation = rad2deg.(eulerangles_from_rotationmatrix(rotation_matrix))
    camera1 = CameraProperties(
            position = position,
            euler_rotation = euler_rotation,
            quaternion_rotation = quaternion_camera_rotation,
            intrinsic = Matrix(I, 3, 3),
    )
    lines_view_1 = get_view(cylinders, camera1)

    error = 0.0

    for i in 1:3
        for j in 1:2
            line = lines_view_1[i, j, :]
            equation = line' * camera1.intrinsic * camera1.rotation_matrix * points_at_infinity[i, :]
            error += equation
        end
    end

    @test isapprox(error, 0.0; atol=1e-6)
end

@testset "calibrated_fixed_position" begin
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)
    position = [10.0, 10.0, -10.0]
    rotation_matrix = lookat_rotation(position, [0.0, 0.0, 0.0])
    quaternion_camera_rotation = QuatRotation(rotation_matrix)
    euler_rotation = rad2deg.(eulerangles_from_rotationmatrix(rotation_matrix))
    camera1 = CameraProperties(
            position = position,
            euler_rotation = euler_rotation,
            quaternion_rotation = quaternion_camera_rotation,
            intrinsic = Matrix(I, 3, 3),
    )
    lines_view_1 = get_view(cylinders, camera1)

    error = 0.0

    for i in 1:3
        for j in 1:2
            line = lines_view_1[i, j, :]
            equation = line' * camera1.intrinsic * camera1.rotation_matrix * points_at_infinity[i, :]
            error += equation
        end
    end

    @test isapprox(error, 0.0; atol=1e-6)
end

@testset "calibrated_random_position" begin
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)
    position, rotation_matrix = random_camera_lookingat_center()
    display(position)
    quaternion_camera_rotation = QuatRotation(rotation_matrix)
    euler_rotation = rad2deg.(eulerangles_from_rotationmatrix(rotation_matrix))
    display(euler_rotation)
    camera1 = CameraProperties(
            position = position,
            euler_rotation = euler_rotation,
            quaternion_rotation = quaternion_camera_rotation,
            intrinsic = Matrix(I, 3, 3),
    )
    lines_view_1 = get_view(cylinders, camera1)

    error = 0.0

    for i in 1:3
        for j in 1:2
            line = lines_view_1[i, j, :]
            equation = line' * camera1.intrinsic * camera1.rotation_matrix * points_at_infinity[i, :]
            error += equation
        end
    end

    @test isapprox(error, 0.0; atol=1e-6)
end

@testset "uncalibrated" begin
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)
    position, rotation_matrix = random_camera_lookingat_center()
    quaternion_camera_rotation = QuatRotation(rotation_matrix)
    euler_rotation = rad2deg.(eulerangles_from_rotationmatrix(rotation_matrix))
    camera1 = CameraProperties(
            position = position,
            euler_rotation = euler_rotation,
            quaternion_rotation = quaternion_camera_rotation,
            intrinsic = intrinsic,
    )
    lines_view_1 = get_view(cylinders, camera1)

    error = 0.0

    for i in 1:3
        for j in 1:2
            line = lines_view_1[i, j, :]
            equation = line' * camera1.intrinsic * camera1.rotation_matrix * points_at_infinity[i, :]
            error += equation
        end
    end

    @test isapprox(error, 0.0; atol=1e-6)
end

@testset "using_camera_matrix" begin
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)
    camera1 = read_camera("../assets/test_scenes/axis_rig/scene.json"; object_path="0.camera")
    lines_view_1 = read_axis_rig_lines("../assets/test_scenes/axis_rig/scene.json"; object_path="0.contours")

    error = 0.0

    for i in 1:3
        for j in 1:2
            line = lines_view_1[i][j]
            equation = line' * camera1.matrix[1:3, 1:3] * points_at_infinity[i, :]
            error += equation
        end
    end

    @test isapprox(error, 0.0; atol=1e-6)
end

@testset "external_uncalibrated" begin
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)
    camera1 = read_camera("../assets/test_scenes/axis_rig/scene.json"; object_path="0.camera")
    lines_view_1 = read_axis_rig_lines("../assets/test_scenes/axis_rig/scene.json"; object_path="0.contours")

    error = 0.0

    for i in 1:3
        for j in 1:2
            line = lines_view_1[i][j]
            equation = line' * camera1.intrinsic * camera1.rotation_matrix * points_at_infinity[i, :]
            error += equation
        end
    end

    @test isapprox(error, 0.0; atol=1e-6)
end

@testset "cayley_rotation" begin
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)
    position, rotation_matrix = random_camera_lookingat_center()
    camera1 = CameraProperties(
            position = position,
            euler_rotation = rad2deg.(eulerangles_from_rotationmatrix(rotation_matrix)),
            quaternion_rotation = QuatRotation(rotation_matrix),
            intrinsic = intrinsic,
    )
    lines_view_1 = get_view(cylinders, camera1)

    quat = Rotations.params(camera1.rotation_matrix)
    quat = quat ./ quat[1]
    quat = quat[2:end]

    rotation = build_rotation_matrix(quat..., false)

    error = 0.0

    for i in 1:3
        for j in 1:2
            line = lines_view_1[i, j, :]
            equation = line' * camera1.intrinsic * rotation * points_at_infinity[i, :]
            error += equation
        end
    end

    @test isapprox(error, 0.0; atol=1e-6)
end

@testset "external_cayley_rotation" begin
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)
    camera1 = read_camera("../assets/test_scenes/axis_rig/scene.json"; object_path="0.camera")
    lines_view_1 = read_axis_rig_lines("../assets/test_scenes/axis_rig/scene.json"; object_path="0.contours")

    quat = Rotations.params(camera1.rotation_matrix)
    quat = quat ./ quat[1]
    quat = quat[2:end]

    rotation = build_rotation_matrix(quat..., false)

    error = 0.0

    for i in 1:3
        for j in 1:2
            line = lines_view_1[i][j]
            equation = line' * camera1.intrinsic * rotation * points_at_infinity[i, :]
            error += equation
        end
    end

    @test isapprox(error, 0.0; atol=1e-6)
end