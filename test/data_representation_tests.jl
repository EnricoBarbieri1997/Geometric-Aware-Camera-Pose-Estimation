using CylindersBasedCameraResectioning.Space: build_rotation_matrix
using CylindersBasedCameraResectioning.Geometry: get_view
using CylindersBasedCameraResectioning.Cylinder: CalibrationRigs, points_at_infinity_dualquadrics
using CylindersBasedCameraResectioning.Camera: CameraProperties, lookat_quaternion, random_camera_lookingat_center
using CylindersBasedCameraResectioning.IO: read_axis_rig_lines, read_camera
using CylindersBasedCameraResectioning.Utils: ≃, eulerangles_from_rotationmatrix

using LinearAlgebra: I
using Rotations
using Random

intrinsics = [
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
        position=position,
        euler_rotation=euler_rotation,
        quaternion_rotation=quaternion_camera_rotation,
        intrinsic=Matrix(I, 3, 3),
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

    @test error ≃ 0.0
end

@testset "calibrated_fixed_position" begin
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)
    position = [10.0, 10.0, -10.0]
    rotation_matrix = lookat_quaternion(position, [0.0, 0.0, 0.0])
    euler_rotation = rad2deg.(eulerangles_from_rotationmatrix(rotation_matrix))
    camera1 = CameraProperties(
        position=position,
        euler_rotation=euler_rotation,
        quaternion_rotation=rotation_matrix,
        intrinsic=Matrix(I, 3, 3),
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

    @test error ≃ 0.0
end

@testset "calibrated_random_position" begin
    Random.seed!(42)
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)
    position, quaternion_camera_rotation = random_camera_lookingat_center()
    euler_rotation = rad2deg.(eulerangles_from_rotationmatrix(quaternion_camera_rotation))
    camera1 = CameraProperties(
        position=position,
        euler_rotation=euler_rotation,
        quaternion_rotation=quaternion_camera_rotation,
        intrinsic=Matrix(I, 3, 3),
    )
    lines_view_1 = get_view(cylinders, camera1)

    error = 0.0

    for i in 1:3
        point_at_infinity_index = (i-1) * 2 + 1
        for j in 1:2
            line = lines_view_1[i, j, :]
            equation = line' * camera1.rotation_matrix * points_at_infinity[point_at_infinity_index, :]
            error += equation
        end
    end

    @test error ≃ 0.0
end

@testset "uncalibrated" begin
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)
    position, rotation_matrix = random_camera_lookingat_center()
    quaternion_camera_rotation = QuatRotation(rotation_matrix)
    euler_rotation = rad2deg.(eulerangles_from_rotationmatrix(rotation_matrix))
    camera1 = CameraProperties(
        position=position,
        euler_rotation=euler_rotation,
        quaternion_rotation=quaternion_camera_rotation,
        intrinsic=intrinsics,
    )
    lines_view_1 = get_view(cylinders, camera1)

    error = 0.0

    for i in 1:3
        point_at_infinity_index = (i-1) * 2 + 1
        for j in 1:2
            line = lines_view_1[i, j, :]
            equation = line' * camera1.intrinsic * camera1.rotation_matrix * points_at_infinity[point_at_infinity_index, :]
            error += equation
        end
    end

    @test error ≃ 0.0
end

@testset "cayley_rotation" begin
    cylinders = CalibrationRigs.axis_rig()
    points_at_infinity, dualquadrics = points_at_infinity_dualquadrics(cylinders)
    position, rotation_matrix = random_camera_lookingat_center()
    camera1 = CameraProperties(
        position=position,
        euler_rotation=rad2deg.(eulerangles_from_rotationmatrix(rotation_matrix)),
        quaternion_rotation=QuatRotation(rotation_matrix),
        intrinsic=intrinsics,
    )
    lines_view_1 = get_view(cylinders, camera1)

    quat = Rotations.params(QuatRotation(camera1.rotation_matrix))
    quat = quat ./ quat[1]
    quat = quat[2:end]

    rotation = build_rotation_matrix(quat..., false)

    error = 0.0

    for i in 1:3
        point_at_infinity_index = (i-1) * 2 + 1
        for j in 1:2
            line = lines_view_1[i, j, :]
            equation = line' * camera1.intrinsic * rotation * points_at_infinity[point_at_infinity_index, :]
            error += equation
        end
    end

    @test error ≃ 0.0
end