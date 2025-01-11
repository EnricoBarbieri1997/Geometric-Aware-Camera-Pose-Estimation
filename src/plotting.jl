module Plotting
    export initfigure, add_2d_axis!, plot_2dpoints, plot_line_2d, Plot3dCameraInput, plot_3dcamera, Plot3dCylindersInput, plot_cylinders_contours, plot_3dcylinders, plot_2dcylinders

    using ..Geometry: Line

    using LinearAlgebra: deg2rad
    using Rotations
    using GLMakie, GLMakie.FileIO

    # f, ax3, ax2

    colors = [:red, :green, :blue, :yellow, :purple, :orange, :pink, :brown]
    linestyles = Dict(
        "dash" => :dash,
        "solid" => :solid,
    )

    function add_2d_axis!()
        index = length(ax2_array) + 1
        row = ceil(Int, index / 3)
        col = index % 3
        if col == 0
            col = 3
        end
        ax = Axis(grid_2d[row, col], autolimitaspect = 1)
        push!(ax2_array, ax)
    end

    function initfigure()
        global f, ax3, grid_2d, ax2_array
        f = Figure(size=(1200, 800))
        ax3 = Axis3(f[1, 1], title = "Cylinders", aspect = :data, perspectiveness = 1.0)
        grid_2d = f[1, 2] = GridLayout()
        Label(grid_2d[:, :, Top()], "Conics")
        ax2_array = []
        add_2d_axis!()
        return f
    end

    struct Plot3dCameraInput
        cameraRotation::Vector{<:Number}
        cameraTranslation::Vector{<:Number}
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
        rotate!(cameraMesh, cameraRotationAxis, cameraRotationAngle)
        translate!(cameraMesh,
            (
                info.cameraTranslation[1],
                info.cameraTranslation[2],
                info.cameraTranslation[3],
            )
        )
    end

    struct Plot3dCylindersInput
        transforms::Vector{Matrix{Float64}}
        radiuses::Vector{Vector{Float64}}
        numberOfCylinders::Int
        cameraProjectionMatrix::Union{Matrix{<:Number}, UndefInitializer}

        function Plot3dCylindersInput(transforms::Vector{Matrix{Float64}}, radiuses::Vector{Vector{Float64}}, numberOfCylinders::Int, cameraProjectionMatrix::Union{Matrix{<:Number}, UndefInitializer} = undef)
            new(transforms, radiuses, numberOfCylinders, cameraProjectionMatrix)
        end
    end
    function plot_3dcylinders(cylindersInfo::Plot3dCylindersInput; axindex = 1)    
        heightLevels = 100
        angles = 100

        z, θ = LinRange(-20, 20, heightLevels), LinRange(0, 2π, angles)
        x = cos.(θ)
        y = sin.(θ)

        for i in 1:cylindersInfo.numberOfCylinders
            radius = cylindersInfo.radiuses[i]
            X = radius[1] * x
            Y = radius[2] * y

            canonicPoints = Array{Float64}(undef, 0, 4)

            for j in 1:heightLevels
                canonicPoints = vcat(canonicPoints, [X Y (z[j] * ones(angles)) ones(angles)])
            end
            points = transpose(cylindersInfo.transforms[i] * canonicPoints')
            points = points ./ points[:, 4]

            lines!(ax3, points[:, 1], points[:, 2], points[:, 3], color = colors[i])

            if (cylindersInfo.cameraProjectionMatrix != undef)
                points2d = [cylindersInfo.cameraProjectionMatrix * point for point in eachrow(points)]
                points2d = [(point ./ point[3]) for point in points2d]
                points2d = hcat(points2d...)'

                lines!(ax2_array[axindex], points2d[:, 1], -points2d[:, 2], color = colors[i])
            end
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
        xs = -50:1:50
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
                xs = -50:1:50
                ys1 = y1.(xs)
                lines!(ax2_array[axindex], xs, -ys1, color = (colors[i], alpha), linestyle=linestyle)
            end
        end
    end
end