using ..Space: RotRad, position_rotation, transformation as create_transform_matrix
using ..Cylinder: CylinderProperties
using ..Geometry: Line, Plane, plane_basis
using ..Camera: CameraProperties

using Reexport
using LinearAlgebra: cross, deg2rad, dot, normalize
using Rotations
using GLMakie
using GLMakie.FileIO
using GeometryBasics

function add_2d_axis!()
    index = length(ax2_array) + 1
    col = ceil(Int, index / 2)
    row = index % 2
    if row == 0
        row = 2
    end
    ax = Axis(grid_2d[row, col], autolimitaspect = 1, aspect = DataAspect(), title="View $index")
    ax.limits[] = ((0, 1080), (-1920, 0))
    push!(ax2_array, ax)
end

function get_or_add_2d_axis!(index)
    if index > length(ax2_array)
        add_2d_axis!()
    end
    return ax2_array[index]
end

function add_slider!(;
    start = 0.0,
    stop = 1.0,
    step = 0.01,
)
    global f
    Box(f[3, :]; height=50)
    slider = Slider(f[3, :]; startvalue=start, range=start:step:stop)
    return slider
end

function add_camera_rotation_axis!()
    index = length(camera_rotation_axes) + 1
    row = ceil(Int, index / 3)
    col = index % 3
    if col == 0
        col = 3
    end
    ax = Axis3(camera_roation_layout[row, col], aspect = :data, perspectiveness = 1.0)
    push!(camera_rotation_axes, ax)
end

function get_or_add_camera_rotation_axis!(index)
    if index > length(camera_rotation_axes)
        add_camera_rotation_axis!()
    end
    return camera_rotation_axes[index]
end

function clean_plots!()
    global f, ax3, ax2_array
    for ax2 in ax2_array
        empty!(ax2)
    end
    for ax_camera in camera_rotation_axes
        empty!(ax_camera)
    end
    empty!(ax3)
end

function initfigure()
    global f, scene3D, ax3, grid_2d, ax2_array, camera_roation_layout, camera_rotation_axes, cameras
    # set_theme!(scale_plot = true)
    f = Figure(size=(1200, 800))
    scene3D = LScene(f[1, 1], scenekw = (camera=cam3d!, show_axis=true))
    ax3 = scene3D.scene
    scatter!(ax3, (0, 0, 0), color = :black, markersize = 10)
    scatter!(ax3, (30, 0, 0), color = :red, markersize = 10)
    scatter!(ax3, (-30, 0, 0), color = :red, markersize = 10, alpha=0.5)
    scatter!(ax3, (0, 30, 0), color = :green, markersize = 10)
    scatter!(ax3, (0, -30, 0), color = :green, markersize = 10, alpha=0.5)
    scatter!(ax3, (0, 0, 30), color = :blue, markersize = 10)
    scatter!(ax3, (0, 0, -30), color = :blue, markersize = 10, alpha=0.5)
    rowsize!(f.layout, 1, Relative(2/3))
    colsize!(f.layout, 1, Relative(2/3))
    grid_2d = f[1, 3] = GridLayout()
    Label(grid_2d[:, :, Top()], "Conics")
    ax2_array = []
    add_2d_axis!()
    camera_roation_layout = f[2, :] = GridLayout()
    camera_rotation_axes = []
    cameras = []

    return f
end

function transform_mesh(mesh::GeometryBasics.Mesh, T::Matrix)
    verts = GeometryBasics.coordinates(mesh)
    new_vertices = [Point3f(T * Vec4f(v..., 1.0)) for v in verts]
    faces = collect(GeometryBasics.faces(mesh))
    return GeometryBasics.Mesh(new_vertices, faces)
end

function plot_frustum!(scene, T; fov_deg=45.0, aspect=1.5, near=0.2, far=10, color=:gray)
    fov = deg2rad(fov_deg)
    h_near = 2 * tan(fov/2) * near
    w_near = h_near * aspect
    h_far = 2 * tan(fov/2) * far
    w_far = h_far * aspect

    # Points in camera space
    points = [
        Vec4f(0, 0, 0, 1),  # camera origin
        Vec4f(-w_near/2, -h_near/2, near, 1),
        Vec4f(w_near/2, -h_near/2, near, 1),
        Vec4f(w_near/2, h_near/2, near, 1),
        Vec4f(-w_near/2, h_near/2, near, 1),
        Vec4f(-w_far/2, -h_far/2, far, 1),
        Vec4f(w_far/2, -h_far/2, far, 1),
        Vec4f(w_far/2, h_far/2, far, 1),
        Vec4f(-w_far/2, h_far/2, far, 1)
    ]

    # Transform to world coordinates
    world_pts = [Point3f(T * p) for p in points]

    # Lines from origin to corners
    for i in 2:5
        lines!(scene, [world_pts[1], world_pts[i]], color=color)
    end

    for i in 6:9
        lines!(scene, [world_pts[1], world_pts[i]], color=color, linestyle=:dot)
    end

    # Edges of near and far planes
    near_ids = [2, 3, 4, 5, 2]
    far_ids = [6, 7, 8, 9, 6]

    lines!(scene, world_pts[near_ids], color=color)
    lines!(scene, world_pts[far_ids], color=color)
end

function plot_3dcamera(camera::CameraProperties, color = :black)
    loaded_mesh = load("./assets/camera.stl")

    scale = 2
    origin = camera.position
    R = Matrix{Float64}(camera.quaternion_rotation)
    T = Matrix(create_transform_matrix(origin, camera.euler_rotation))

    # Plot axes
    for (i, color) in enumerate((:red, :green, :blue))
        endpoint = origin .+ scale * R[:, i]
        lines!(ax3, [origin[1], endpoint[1]], [origin[2], endpoint[2]], [origin[3], endpoint[3]],
               color=color, linewidth=2)
    end

    fx = camera.intrinsic[1, 1]
    fov_x_rad = 2 * atan(1080 / (2 * fx))
    fov_x_deg = rad2deg(fov_x_rad)
    aspect = camera.intrinsic[2, 2] / camera.intrinsic[1, 1]
    plot_frustum!(ax3, T; fov_deg = fov_x_deg, aspect = aspect, color = color)

    mesh_transformed = transform_mesh(loaded_mesh, T)
    mesh!(ax3, mesh_transformed, color=:gray)
    push!(cameras, Dict(
        :camera => camera,
        :mesh => mesh_transformed,
    ))
    index = length(cameras)
    text!(ax3, "Camera $(index)", position=Point3f(camera.position), align = (:center, :center), color=color)
end

function plot_3dcamera_rotation(camera::CameraProperties; color = :black, axindex = nothing)
    ax = if !isnothing(axindex) camera_rotation_axes[axindex] else ax3 end
    cameraModel = load("./assets/camera.stl")
    cameraMesh = mesh!(
        ax,
        cameraModel,
        color = color,
    )
    cameraRotationRad = deg2rad.(camera.euler_rotation)
    cameraRotation = RotRad(cameraRotationRad...)
    cameraRotationAxis = rotation_axis(cameraRotation)
    cameraRotationAngle = rotation_angle(cameraRotation)
    rotate!(cameraMesh, cameraRotationAxis, cameraRotationAngle)
end

function plot_plane!(π::Vector{<:Real}; corners=nothing, origin=nothing, color=:gray)
    # Step 1: Transform the plane
    πt = π
    n = πt[1:3]
    d = πt[4]
    n = normalize(n)

    # Step 2: Get a point on the plane
    p0 = something(origin, -d * n)[1:3]

    u, v = plane_basis(Plane(p0, n))
    u = u[1:3]
    v = v[1:3]

    size = 4.0

    # Step 4: Generate 4 points (quad corners)
    corners = something(corners, Point3f[
        p0 + size*u + size*v,
        p0 - size*u + size*v,
        p0 - size*u - size*v,
        p0 + size*u - size*v,
    ])

    # Step 5: Define faces as a quad (two triangles)
    faces = [
        TriangleFace(1, 2, 3),
        TriangleFace(1, 3, 4),
        TriangleFace(1, 3, 2),
        TriangleFace(1, 4, 3),
    ]

    # Step 6: Create mesh and render with double-sided shading
    mesh!(ax3, GeometryBasics.Mesh(corners, faces),
          color = color, transparency = true,
          shading = NoShading)
end

function plot_3dcylinders(cylinders::Vector{CylinderProperties}; axindex = 1)
    for (i, cylinder) in enumerate(cylinders)
        P0, _ = position_rotation(cylinder.transform)                 # Base point of the cylinder
        scatter!(ax3, [
            Point3f(P0[1], P0[2], P0[3])
        ]; color = colors[i], markersize = 40)
        v = normalize(cylinder.singular_point[1:3])       # Axis direction (must be normalized)
        r = cylinder.radiuses[1]                              # Cylinder radius
        height = 10.0                         # Total height of the cylinder
        angle_step = 5                      # Angular step in degrees
        height_step = 0.1               # Height step
    
        a = abs(v[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
        u = normalize(cross(v, a))
        w = cross(v, u)
    
        # --- Generate points ---
        start_height = -height + height_step
        end_height = height - height_step
        height_sections = start_height:height_step:end_height
        angles = 0:angle_step:355
        number_of_height_sections = size(height_sections)[1]
        number_of_angle_sections = size(angles)[1]
        points = zeros(Float64, 3, number_of_height_sections * number_of_angle_sections)

        for (h_i, h) in enumerate(height_sections)
            for (deg_i, deg) in enumerate(angles)
                θ = deg2rad(deg)
                point = P0 .+ h .* v .+ r*cos(θ).*u .+ r*sin(θ).*w
                index = (h_i - 1) * number_of_angle_sections + deg_i
                points[:, index] = point
            end
        end

        lines!(ax3, points[1, :], points[2, :], points[3, :]; color = colors[i])

        plane = inv(cylinder.transform') * [0, 0, 1, 0]
        plot_plane!(plane; origin=P0, color=colors[i])
    end
end

function plot_2dpoints(points; axindex = 1)
    for (i, point) in enumerate(points)
        if length(point) == 3
            point = point ./ point[3]
        end
        scatter!(ax2_array[axindex], (point[1], -point[2]), color = colors[i])
    end
end

function plot_line_2d(line:: Line; color = :black, linestyle = :solid, axindex = 1)
    slope = line.direction[2] / line.direction[1]
    intercept = line.origin[2] - slope * line.origin[1]

    y = function (x) return slope * x + intercept end
    xs = 0:1:1080
    ys = y.(xs)
    lines!(ax2_array[axindex], xs, ys, color = color, linestyle=linestyle)
end

function plot_cylinders_contours(contours::Vector{Vector{Line}}; linestyle = :solid)
    for (i, contour) in enumerate(contours)
        for line in contour
            plot_line_2d(line, color = colors[i], linestyle = linestyle)
        end
    end
end

function plot_2dcylinders(conic_contours; linestyle = :solid, alpha = 1, axindex = 1)
    y = function (x, l) return (-(l[1] * x + l[3]) / l[2]) end
    for i in 1:(size(conic_contours)[1])
        for j in 1:(size(conic_contours)[2])
            line = conic_contours[i, j, :]
            if line == [0, 0, 0] continue end
            y1 = function (x) return y(x, line) end
            xs = 0:1:1080
            ys1 = y1.(xs)
            lines!(ax2_array[axindex], xs, -ys1, color = (colors[i], alpha), linestyle=linestyle)
        end
    end
end

function plot_image_background(img; axindex = 1)
    _, h = size(img)
    image!(ax2_array[axindex], img;
        transformation = (scale = Vec3f(1), translation = Vec3f(0, -1 * h, -1)),
    )
end
