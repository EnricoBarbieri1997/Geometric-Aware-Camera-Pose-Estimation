using CylindersBasedCameraResectioning.Geometry
using CylindersBasedCameraResectioning: sample_camera
using CylindersBasedCameraResectioning.Cylinder: CylinderProperties, CalibrationRigs
using CylindersBasedCameraResectioning.Camera: CameraProperties, random_camera_lookingat_center
using CylindersBasedCameraResectioning.IO: read_axis_rig_lines, read_camera
using CylindersBasedCameraResectioning.Utils: rand_in_range

using LinearAlgebra: norm, normalize
using Random

@testset "issame_line" begin
    line₁ = Line([0, 0, 0], [1, 1, 1])
    line₂ = Line([1, 1, 1], [2, 2, 2])
    @test issame_line(line₁, line₂)
end

@testset "project_point_into_line" begin
    point = [1, 1, 1]
    line = Line([0, 0, 0], [0, 0, 1])
    projected_point = project_point_into_line(point, line)
    @test projected_point ≈ [0, 0, 1]
end

@testset "project_point_into_plane" begin
    point = [1, 1, 1]
    plane = Plane([0, 0, 0], [0, 0, 1])
    projected_point = project_point_into_plane(point, plane)
    @test projected_point ≈ [1, 1, 0]
end

@testset "get_tangentpoints_circle_point" begin
    circle = Circle([1, 3, 2], 4, [0, -2, -1])
    point = [-6, 3, 2]
    correctpoint₁ = [-1.29, 4.47, -0.94]
    correctpoint₂ = [-1.29, 1.53, 4.94]
    tangentpoint₁, tangentpoint₂ = get_tangentpoints_circle_point(circle, point)
    tangentpoint₁ = round.(tangentpoint₁; digits=2)
    tangentpoint₂ = round.(tangentpoint₂; digits=2)

    @test (tangentpoint₁ == correctpoint₁ && tangentpoint₂ == correctpoint₂) ||
          (tangentpoint₂ == correctpoint₁ && tangentpoint₁ == correctpoint₂)
end

@testset "get_cylinder_contours" begin
    cylinders = CalibrationRigs.axis_rig()
    camera = read_camera("../assets/test_scenes/axis_rig/scene.json"; object_path="0.camera")
    target_lines = read_axis_rig_lines("../assets/test_scenes/axis_rig/scene.json"; object_path="0.contours")

    for (i, cylinder) in enumerate(cylinders)
        lines = get_cylinder_contours(
            cylinder,
            camera
        )
        for (j, _line) in enumerate(lines)
            line = _line ./ _line[3]
            target_line = target_lines[i][j]
            target_line = target_line ./ target_line[3]

            @test isapprox(norm(line - target_line), 0.0; atol=1e-2)
        end
    end
end

@testset "compute_Hinf" begin
    Random.seed!(84564)

    # 4 random vanishing points (3D directions, homogeneous with w=0)
    n_points = 4
    vanishing_points = hcat([[normalize(randn(3)); 0.0] for _ in 1:n_points]...)

    # Shared intrinsics for both cameras
    intrinsics = [
        rand_in_range(2500.0, 2700.0) 0.0 rand_in_range(950.0, 970.0);
        0.0 rand_in_range(1400.0, 1600.0) rand_in_range(530.0, 550.0);
        0.0 0.0 1.0
    ]

    # Create 2 random cameras
    cameras = CameraProperties[]
    for _ in 1:2
        position, rotation = random_camera_lookingat_center()
        camera = CameraProperties()
        camera.position = position
        camera.quaternion_rotation = rotation
        camera.intrinsic = intrinsics
        push!(cameras, camera)
    end

    # Project vanishing points through each camera
    function project_vps(camera, vps)
        n = size(vps, 2)
        pts_2d = zeros(n, 2)
        for i in 1:n
            projected = camera.matrix * vps[:, i]
            pts_2d[i, :] = projected[1:2] ./ projected[3]
        end
        return pts_2d
    end

    view1 = project_vps(cameras[1], vanishing_points)
    view2 = project_vps(cameras[2], vanishing_points)

    # Ground truth H_∞ = K₂ * R₂ * R₁' * K₁⁻¹
    K1, K2 = cameras[1].intrinsic, cameras[2].intrinsic
    R1, R2 = cameras[1].rotation_matrix, cameras[2].rotation_matrix
    Hinf_true = K2 * R2 * R1' * inv(K1)
    Hinf_true = Hinf_true ./ Hinf_true[3, 3]

    # Compute H_∞ using DLT
    Hinf_dlt = compute_Hinf(view1, view2)

    # Verify both have low reprojection errors
    for i in 1:n_points
        p1_h = [view1[i, :]; 1.0]

        # Ground truth error
        p2_pred_true = Hinf_true * p1_h
        p2_pred_true_inhom = p2_pred_true[1:2] ./ p2_pred_true[3]
        error_true = norm(p2_pred_true_inhom - view2[i, :])
        @test error_true < 1e-10

        # DLT error
        p2_pred_dlt = Hinf_dlt * p1_h
        p2_pred_dlt_inhom = p2_pred_dlt[1:2] ./ p2_pred_dlt[3]
        error_dlt = norm(p2_pred_dlt_inhom - view2[i, :])
        @test error_dlt < 1e-10
    end
end