using CylindersBasedCameraResectioning.Cylinder: CylinderProperties, cylinder_from_center_axis_radius
using CylindersBasedCameraResectioning.IO: read_cylinders_yup_to_zup

using LinearAlgebra: dot, norm
using JSON

function zup_to_yup(vec::Vector{<:Number})
    x, y, z = vec
    return Float64[x, z, -y]
end

@testset "read_cylinders_yup_to_zup -> CylinderProperties" begin
    gt_cylinder = cylinder_from_center_axis_radius([1.2, -3.4, 5.6], [0.3, -0.4, 0.5], 2.1)
    gt_center = gt_cylinder.transform[1:3, 4]
    gt_axis = gt_cylinder.singular_point[1:3]
    gt_axis = gt_axis / norm(gt_axis)
    gt_radius = gt_cylinder.radiuses[1]

    cylinder_json = Dict(
        "cylinders" => [Dict(
            "center" => zup_to_yup(collect(gt_center)),
            "axis" => zup_to_yup(collect(gt_axis)),
            "radius" => gt_radius,
        )],
    )

    filepath = joinpath(@__DIR__, "data", "cylinders_y_up.json")
    mkpath(dirname(filepath))
    open(filepath, "w") do io
        JSON.print(io, cylinder_json)
    end

    parsed = read_cylinders_yup_to_zup(filepath)
    @test length(parsed) == 1

    parsed_cylinder = parsed[1]
    reconstructed = cylinder_from_center_axis_radius(parsed_cylinder.center, parsed_cylinder.axis, parsed_cylinder.radius)

    reconstructed_axis = reconstructed.singular_point[1:3]
    reconstructed_axis = reconstructed_axis / norm(reconstructed_axis)

    @test isapprox(reconstructed.transform[1:3, 4], gt_center; atol=1e-10)
    @test isapprox(reconstructed.radiuses[1], gt_radius; atol=1e-10)
    @test isapprox(abs(dot(reconstructed_axis, gt_axis)), 1.0; atol=1e-10)
    @test isapprox(reconstructed.matrix, gt_cylinder.matrix; atol=1e-8)
    @test isapprox(reconstructed.dual_matrix, gt_cylinder.dual_matrix; atol=1e-8)
end
