module Minimality

using LinearAlgebra: cross, dot, norm, normalize, diagm, inv, hcat, vcat
using ..Geometry: Circle, Line, get_cylinder_contours
using ..Space: build_rotation_matrix
using ..Camera: CameraProperties, build_camera_matrix
using ..Cylinder: CylinderProperties, standard_and_dual as standard_and_dual_cylinder
using ..Utils: eulerangles_from_rotationmatrix

# using GLMakie: lines! as draw_lines!, Figure, Axis, DataAspect, xlims!, ylims!
using ..Plotting
using HomotopyContinuation

function toZero(val)
    if val < 1E-6
        return 0
    end
    return val
end

# function get_tangentpoints_circle_point(circle::Circle, point::Vector{<:Number})
#     variation = point - circle.center
#     d = norm(variation)
#     r = circle.radius

#     if d <= r
#         throw(TangentLineNotFound("The point is inside the circle"))
#     end
#     if (d / r - 1) < 1E-8
#         throw(TangentLineNotFound("The point is on the circle"))
#     end

#     R = sqrt(d^2 - r^2)
#     rho = r / d
#     ad = rho^2 * d
#     bd = rho * R / d * d
#     axis = normalize(circle.axis)

#     adv = normalize(variation)
#     if !isnothing(axis)
#         bdv = normalize(cross(axis, adv))
#     else
#         bdv = normalize(adv .* [1, -1])
#     end

#     T1 = circle.center + ad * adv + bd * bdv
#     T2 = circle.center + ad * adv - bd * bdv

#     # T1 = rationalize.(toZero.(T1))
#     # T2 = rationalize.(toZero.(T2))
#     return T1, T2

#     throw(TangentLineNotFound("No tangent line possible"))
# end

# function build_camera_matrix(intrinsic, rotation, translation; use_rotation_as_is=false, use_translation_as_is=false)
#     r₁ = rotation
#     if (!use_rotation_as_is)
#         r₁ = r₁' # inv(r)
#     end
#     t₁ = translation
#     if (!use_translation_as_is)
#         t₁ = -r₁ * translation
#     end
#     m = hcat(intrinsic, zeros(Rational, 3)) * vcat(hcat(r₁, t₁), [0 0 0 1])
#     return m ./ m[3, 4]
# end

# function project_point_into_line(point::Vector{<:Number}, line::Line)::Vector{<:Number}
#     # direction = rationalize.(line.direction / norm(line.direction))
#     direction = line.direction / norm(line.direction)
#     origin = line.origin
#     v = point - origin
#     return origin + dot(v, direction) * direction
# end

# function get_cylinder_contours(cylinder_axis, camera_position, camera_matrix)
#     radius = 1
#     cylinder_axis = cylinder_axis[1:3]
#     circlecenter = project_point_into_line(camera_position, Line([0, 0, 0], cylinder_axis))
#     tangentpoint₁, tangentpoint₂ = get_tangentpoints_circle_point(
#         Circle(circlecenter, radius, cylinder_axis),
#         camera_position
#     )

#     projected_tangentpoint₁ = (camera_matrix * [tangentpoint₁; 1])
#     projected_tangentpoint₂ = (camera_matrix * [tangentpoint₂; 1])
#     projected_cylinderaxis = (camera_matrix * [cylinder_axis; 0])

#     contour₁ = cross(projected_tangentpoint₁, projected_tangentpoint₁ + projected_cylinderaxis)
#     contour₂ = cross(projected_tangentpoint₂, projected_tangentpoint₂ + projected_cylinderaxis)

#     return (contour₁, contour₂)
# end

pi = 22 // 7

function equations()
    real_intrinsics::Matrix{Rational} = [
        1500 0 540;
        0 2667 960;
        0 0 1
    ]
    cylinders_transform_matrices = [
        [
            1 0 0 0;
            0 1 0 0;
            0 0 1 0;
            0 0 0 1
        ],
        [
            1 0 0 0;
            0 0 1 0;
            0 -1 0 0;
            0 0 0 1
        ],
        [
            0 0 1 0;
            0 1 0 0;
            -1 0 0 0;
            0 0 0 1
        ]
    ]

    cylinders_dual_matrices = []
    for i in 1:3
        transform_matrix = cylinders_transform_matrices[i]
		canonical_cylinder = diagm([1, 1, 0, -1])
        dual_canonical_cylinder = zeros(4, 4)
		dual_canonical_cylinder[[1, 2, 4], [1, 2, 4]] .= inv(canonical_cylinder[[1, 2, 4], [1, 2, 4]])
        dual_cylinder_matrix = transform_matrix * dual_canonical_cylinder * transform_matrix'
        push!(cylinders_dual_matrices, dual_cylinder_matrix)
    end

    cameras_positions = [
        [10, 10, 10],
        [-10, 10, 10]
    ]
    # camera_rotations_raw = [
    #     [271 // 12, 123 // 8, -11 // 16], # [3π / 4, 0, 5π / 4],
    #     [13//25, 44//35, -70//29]
    # ]
    camera_rotations_raw = [
        [1, 1, 1],
        [1, 3, 1]
    ]
    camera_rotations = [
        build_rotation_matrix(camera_rotations_raw[1]..., true),
        build_rotation_matrix(camera_rotations_raw[2]..., true),
    ]
    display("Camera rotations")
    display(camera_rotations)
    display("-----------------------")
    camera_matrices = [
        build_camera_matrix(real_intrinsics, camera_rotations[1], cameras_positions[1]),
        build_camera_matrix(real_intrinsics, camera_rotations[2], cameras_positions[2])
    ]
    display("Camera matrices")
    display(camera_matrices)
    display("-----------------------")
    points_at_infinity = []
    for cylinder_index in 1:3
        push!(points_at_infinity, cylinders_transform_matrices[cylinder_index] * [0; 0; 1; 0])
    end
    display("Points at infinity")
    display(points_at_infinity)
    display("-----------------------")
    # lines = zeros(Rational{Int}, 2, 3, 2, 4)
    lines = zeros(Float64, 2, 3, 2, 3)

    f = initfigure()
    axx = [
        get_or_add_2d_axis!(1),
        get_or_add_2d_axis!(2)
    ]

    for camera_index in 1:2
        ca = CameraProperties()
        ca.position = cameras_positions[camera_index]
        ca.quaternion_rotation = camera_rotations[camera_index]
        ca.euler_rotation = rad2deg.(eulerangles_from_rotationmatrix(camera_rotations[camera_index]))
        ca.intrinsic = real_intrinsics
        plot_3dcamera(ca)
        cylinders_3d_to_view::Vector{CylinderProperties} = []
        for cylinder_index in 1:3
            cy = CylinderProperties()
            cy.transform = cylinders_transform_matrices[cylinder_index]
            cy.radiuses = [1, 1]
            standard, dual, singularpoint = standard_and_dual_cylinder(cy.transform, cy.radiuses)
            cy.matrix = standard
            cy.dual_matrix = dual
            cy.singular_point = singularpoint
            push!(cylinders_3d_to_view, cy)

            line1, line2 = get_cylinder_contours(
                cy,
                ca
            )
            # line1, line2 = get_cylinder_contours(
            #     points_at_infinity[cylinder_index],
            #     cameras_positions[camera_index],
            #     camera_matrices[camera_index]
            # )
            lines[camera_index, cylinder_index, 1, :] = line1
            lines[camera_index, cylinder_index, 2, :] = line2
        end

        if (camera_index == 1)
            plot_3dcylinders(cylinders_3d_to_view)
        end
    end
    display("Lines")

    for camera_index in 1:2
        contours = zeros(Float64, 3, 2, 3)
        for cylinder_index in 1:3
            for line_index in 1:2
                display("Camera $(camera_index), Cylinder $(cylinder_index), Line $(line_index)")
                line = lines[camera_index, cylinder_index, line_index, :]
                display(line)
                # if cylinder_index == 2
                #     err = (line' * camera_matrices[camera_index] * cylinders_dual_matrices[cylinder_index] * camera_matrices[camera_index]' * line)
                #     @assert err < 1E-1 "(9) Line $(line_index) of conic $(cylinder_index) is tangent to the projected cylinder. $(err)"
                # end
                contours[cylinder_index, line_index, :] .= line
            end
        end
        plot_2dcylinders(contours; axindex=camera_index)
    end
    display(f)
    display("---------------------------")
    system_to_solve = []
    parameters::Vector{HomotopyContinuation.ModelKit.Variable} = []

    @var fx, fy, skew, cx, cy

    intrinsic_variables = [
        fx skew cx;
        0 fy cy;
        0 0 1
    ]

    @var r1x, r1y, r1z, r2x, r2y, r2z

    r1 = build_rotation_matrix(r1x, r1y, r1z, false)
    r2 = build_rotation_matrix(r2x, r2y, r2z, false)
    rotations = [
        r1,
        r2
    ]

    variables = [fx, fy, skew, cx, cy, r1x, r1y, r1z, r2x, r2y, r2z]

    for view_index in 1:2
        for cylinder_index in 1:3
            for line_index in 1:2
                equation = lines[view_index, cylinder_index, line_index, 1:3]' * intrinsic_variables * rotations[view_index] * points_at_infinity[cylinder_index][1:3]
                push!(system_to_solve, equation)
            end
        end
    end

    system_to_solve = system_to_solve[1:11]
    display("System to solve")
    display(system_to_solve)
    display("----------------")
    a = System(system_to_solve; variables)
    display(a)
    display("----------------")

    display("Verification")
    b = evaluate(a, [
        real_intrinsics[1, 1],
        real_intrinsics[2, 2],
        real_intrinsics[1, 2],
        real_intrinsics[1, 3],
        real_intrinsics[2, 3],
        camera_rotations_raw[1][1],
        camera_rotations_raw[1][2],
        camera_rotations_raw[1][3],
        camera_rotations_raw[2][1],
        camera_rotations_raw[2][2],
        camera_rotations_raw[2][3]
    ])
    display(b)
end

end