using CylindersBasedCameraResectioning.Geometry
using CylindersBasedCameraResectioning: sample_camera
using CylindersBasedCameraResectioning.Cylinder: CylinderProperties, CalibrationRigs
using CylindersBasedCameraResectioning.IO: read_axis_rig_lines, read_camera

using LinearAlgebra: norm

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

            @test norm(line - target_line) ≈ 0.0
        end
    end
end