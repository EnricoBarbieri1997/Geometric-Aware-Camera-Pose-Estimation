using ..Space: position_rotation
using ..Cylinder: CylinderProperties
using ..Geometry: Line

using Reexport
using LinearAlgebra: cross, deg2rad, normalize
using Rotations
using GLMakie
using GLMakie.FileIO

function add_2d_axis!()
    index = length(ax2_array) + 1
    row = ceil(Int, index / 2)
    col = index % 2
    if col == 0
        col = 2
    end
    ax = Axis(grid_2d[row, col], autolimitaspect = 1)
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
    global f, ax3, grid_2d, ax2_array, camera_roation_layout, camera_rotation_axes
    f = Figure(size=(1200, 800))
    ax3 = LScene(f[1, 1], scenekw = (show_axis = true,))
    # ax3.limits[] = ((-30, 30), (-30, 30), (-30, 30))
    grid_2d = f[1, 2] = GridLayout()
    Label(grid_2d[:, :, Top()], "Conics")
    ax2_array = []
    add_2d_axis!()
    camera_roation_layout = f[2, :] = GridLayout()
    camera_rotation_axes = []
    return f
end

function plot_3dcamera(info::Plot3dCameraInput, color = :black)
    cameraModel = load("./assets/camera.stl")
    cameraMesh = mesh!(
        ax3,
        cameraModel,
        color = color,
    )
    cameraRotationRad = deg2rad.(info.cameraRotation)
    cameraRotation = RotXYZ(cameraRotationRad...)
    cameraRotationAxis = rotation_axis(cameraRotation)
    cameraRotationAngle = rotation_angle(cameraRotation)
    scale!(cameraMesh, (1/1, 1/1, 1/1))
    rotate!(cameraMesh, cameraRotationAxis, cameraRotationAngle)
    translate!(cameraMesh,
        (
            info.cameraTranslation[1],
            info.cameraTranslation[2],
            info.cameraTranslation[3],
        )
    )
end

function plot_3dcamera_rotation(info::Plot3dCameraInput; color = :black, axindex = nothing)
    ax = if !isnothing(axindex) camera_rotation_axes[axindex] else ax3 end
    cameraModel = load("./assets/camera.stl")
    cameraMesh = mesh!(
        ax,
        cameraModel,
        color = color,
    )
    cameraRotationRad = deg2rad.(info.cameraRotation)
    cameraRotation = RotXYZ(cameraRotationRad...)
    cameraRotationAxis = rotation_axis(cameraRotation)
    cameraRotationAngle = rotation_angle(cameraRotation)
    rotate!(cameraMesh, cameraRotationAxis, cameraRotationAngle)
end

function plot_3dcylinders(cylinders::Vector{CylinderProperties}; axindex = 1)
    for (i, cylinder) in enumerate(cylinders)
        P0, _ = position_rotation(cylinder.matrix)                # Base point of the cylinder
        v = normalize(cylinder.singular_point[1:3])       # Axis direction (must be normalized)
        r = cylinder.radiuses[1]                              # Cylinder radius
        height = 10.0                         # Total height of the cylinder
        angle_step = 5                      # Angular step in degrees
        height_step = 0.1               # Height step
    
        a = abs(v[1]) < 0.9 ? [1.0, 0.0, 0.0] : [0.0, 1.0, 0.0]
        u = normalize(cross(v, a))
        w = cross(v, u)
    
        # --- Generate points ---
        points = zeros(Float64, 3, Int(ceil(360/angle_step) * ceil(2*height/height_step)))

        for h in -height:height_step:(height - height_step)
            for deg in 0:angle_step:355
                θ = deg2rad(deg)
                point = P0 .+ h .* v .+ r*cos(θ).*u .+ r*sin(θ).*w
                index = Int(ceil(deg/angle_step) + 1 + ceil((h+height)/height_step) * ceil(360/angle_step))
                points[:, index] = point
            end
        end

        lines!(ax3, points[1, :], points[2, :], points[3, :]; color = colors[i])
    end
end

function plot_2dpoints(points; axindex = 1)
    for (i, point) in enumerate(points)
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
            lines!(ax2_array[axindex], -xs, -ys1, color = (colors[i], alpha), linestyle=linestyle)
        end
    end
end

function plot_image_background(img; axindex = 1)
    _, h = size(img)
    image!(ax2_array[axindex], img;
        transformation = (scale = Vec3f(1), translation = Vec3f(0, -1 * h, -1)),
    )
end
