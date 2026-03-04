module Cylinder
export CylinderProperties, standard_and_dual, points_at_infinity_dualquadrics, cylinder_from_center_axis_radius

using ..Space, ..Utils

using LinearAlgebra, Random
using Rotations

mutable struct CylinderProperties
    euler_rotation::Vector{Number}
    radiuses::Array{<:Number,1}
    matrix::Matrix{<:Number}
    singular_point::Vector{Number}
    dual_matrix::Matrix{<:Number}
    transform::Matrix{<:Number}

    CylinderProperties() = new()
end

function standard_and_dual(
    transform_matrix::Matrix{<:Number},
    radius::Union{Vector{<:Number},Array{<:Number,1}}=[1, 1]
)
    inverseRadiusSquareX = 1 / (radius[1]^2)
    inverseRadiusSquareY = 1 / (radius[2]^2)
    canonicalCylinder = diagm([inverseRadiusSquareX, inverseRadiusSquareY, 0, -1])

    cylinder = canonicalCylinder
    try
        cylinder = inv(transform_matrix') * canonicalCylinder * inv(transform_matrix)
    catch e
        @warn "Matrix inversion failed in standard_and_dual function. Returning canonical cylinder. Error: $e"
    end

    dualCanonicalCylinderMatrix = zeros(4, 4)
    dualCanonicalCylinderMatrix[[1, 2, 4], [1, 2, 4]] .= inv(canonicalCylinder[[1, 2, 4], [1, 2, 4]])

    dualCylinderMatrix = transform_matrix * dualCanonicalCylinderMatrix * transform_matrix'
    dualCylinderSingularPoint = transform_matrix * [0; 0; 1; 0]

    return cylinder, dualCylinderMatrix, dualCylinderSingularPoint
end

function points_at_infinity_dualquadrics(cylinders)
    l = length(cylinders) * 2
    points_at_infinity::Matrix{Float64} = zeros(Float64, l, 3)
    for (i, cylinder) in enumerate(cylinders)
        doubled_index = (i - 1) * 2 + 1
        points_at_infinity[doubled_index, :] = normalize(cylinder.singular_point[1:3])
        points_at_infinity[doubled_index+1, :] = normalize(cylinder.singular_point[1:3])
    end
    dualquadrics::Array{Float64,3} = zeros(Float64, l, 4, 4)
    for (i, cylinder) in enumerate(cylinders)
        doubled_index = (i - 1) * 2 + 1
        dualquadrics[doubled_index, :, :] = cylinder.dual_matrix ./ cylinder.dual_matrix[4, 4]
        dualquadrics[doubled_index+1, :, :] = cylinder.dual_matrix ./ cylinder.dual_matrix[4, 4]
    end

    return points_at_infinity, dualquadrics
end

function cylinder_from_center_axis_radius(
    center::Vector{<:Number},
    axis::Vector{<:Number},
    radius::Number,
)::CylinderProperties
    axis_direction = Float64.(axis[1:3])
    axis_norm = norm(axis_direction)
    if axis_norm == 0
        throw(ArgumentError("Axis vector must be non-zero."))
    end
    axis_direction = axis_direction / axis_norm

    x_axis = Float64.(get_any_perpendicular(axis_direction))
    y_axis = Float64.(normalize(cross(axis_direction, x_axis)))
    rotation = RotMatrix3(hcat(x_axis, y_axis, axis_direction))
    transform = transformation(Float64.(center[1:3]), rotation)

    radius_float = Float64(radius)
    standard, dual, singularpoint = standard_and_dual(transform, [radius_float, radius_float])

    cylinder = CylinderProperties()
    cylinder.euler_rotation = eulerangles_from_rotationmatrix(rotation)
    cylinder.radiuses = [radius_float, radius_float]
    cylinder.matrix = standard
    cylinder.singular_point = singularpoint
    cylinder.dual_matrix = dual
    cylinder.transform = transform
    return cylinder
end

module CalibrationRigs
using ..Cylinder: CylinderProperties, standard_and_dual
using ....CylindersBasedCameraResectioning: ASSERTS_ENABLED
using ....Space: transformation_euler
using ....Utils: rand_in_range

using LinearAlgebra: normalize
using Random

function cylinder(position, euler_rotation)
    cylinder = CylinderProperties()
    cylinder.euler_rotation = euler_rotation
    cylinder.transform = transformation_euler(position, cylinder.euler_rotation)
    cylinder.radiuses = [1.0, 1.0]

    axis = cylinder.transform * [0; 0; 1; 0]
    axis = axis ./ axis[3]
    axis = axis[1:3]

    standard, dual, singularpoint = standard_and_dual(cylinder.transform, cylinder.radiuses)
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

function origin_centered_cylinder(euler_rotation)
    return cylinder([0.0, 0.0, 0.0], euler_rotation)
end

function axis_rig()
    cylinders = []
    # Cylinder X
    push!(cylinders, origin_centered_cylinder([0.0, 90.0, 0.0]))
    # Cylinder Y
    push!(cylinders, origin_centered_cylinder([90.0, 0.0, 0.0]))
    # Cylinder Z
    push!(cylinders, origin_centered_cylinder([0.0, 0.0, 0.0]))
    return cylinders
end

function arbitrary_rig()
    Random.seed!(2300)
    cylinders = []
    for i in 1:3
        dir = normalize(rand_in_range(-1.0, 1.0, 3))
        position = dir * rand_in_range(0.1, 3.0)
        euler_rotation = rand_in_range(0.0, 360.0, 3)
        push!(cylinders, cylinder(position, euler_rotation))
    end
    return cylinders
end

function arbitrary_rig_four()
    Random.seed!(2300)
    cylinders = []
    for i in 1:4
        dir = normalize(rand_in_range(-1.0, 1.0, 3))
        position = dir * rand_in_range(0.1, 3.0)
        euler_rotation = rand_in_range(0.0, 360.0, 3)
        push!(cylinders, cylinder(position, euler_rotation))
    end
    return cylinders
end
end
end
