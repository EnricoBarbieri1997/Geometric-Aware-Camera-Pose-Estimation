module Plotting
    export initFigure, plot2DPoints, Plot3DCameraInput, plot3DCamera, Plot3DCylindersInput, plot3DCylinders, plot2DCylinders

    using LinearAlgebra: deg2rad
    using Rotations
    using GLMakie, GLMakie.FileIO

    # f, ax3, ax2

    colors = [:red, :green, :blue, :yellow, :purple, :orange, :pink, :brown]

    function initFigure()
        global f, ax3, ax2
        f = Figure(size=(1200, 800))
        ax3 = Axis3(f[1, 1], title = "Cylinders", aspect = :data, perspectiveness = 1.0)
        ax2 = Axis(f[1, 2], title = "Conics", autolimitaspect = 1)
        return f
    end

    struct Plot3DCameraInput
        cameraRotation::Tuple{Float64, Float64, Float64}
        cameraTranslation::Tuple{Float64, Float64, Float64}
    end
    function plot3DCamera(info::Plot3DCameraInput)    
        cameraModel = load("./camera.stl")
        cameraMesh = mesh!(
            ax3,
            cameraModel,
        )
        cameraRotationRad = deg2rad.(info.cameraRotation)
        cameraRotation = RotXYZ(cameraRotationRad...)
        cameraRotationAxis = rotation_axis(cameraRotation)
        cameraRotationAngle = rotation_angle(cameraRotation)
        rotate!(cameraMesh, cameraRotationAxis, cameraRotationAngle)
        translate!(cameraMesh, info.cameraTranslation)
    end

    struct Plot3DCylindersInput
        transforms::Array{Matrix{Float64}}
        radiuses::Array{Tuple{Number, Number}}
        numberOfCylinders::Int
        cameraProjectionMatrix::Union{Matrix{Float64}, UndefInitializer}

        function Plot3DCylindersInput(transforms::Array{Matrix{Float64}}, radiuses::Array{Tuple{Number, Number}}, numberOfCylinders::Int, cameraProjectionMatrix::Union{Matrix{Float64}, UndefInitializer} = undef)
            new(transforms, radiuses, numberOfCylinders, cameraProjectionMatrix)
        end
    end
    function plot3DCylinders(cylindersInfo::Plot3DCylindersInput)    
        heightLevels = 100
        angles = 100
    
        z, θ = LinRange(-20, 20, heightLevels), LinRange(0, 2π, angles)
        x = cos.(θ)
        y = sin.(θ)
    
        for i in 1:cylindersInfo.numberOfCylinders
            radius = cylindersInfo.radiuses[i]
            X = radius[1] * x
            Y = radius[2] * y
    
            canonicPoints = []
    
            for j in 1:heightLevels
                canonicPoints = vcat(canonicPoints, [X Y (z[j] * ones(angles)) ones(angles)])
            end
            points = transpose(cylindersInfo.transforms[i] * canonicPoints')
            points = points ./ points[:, 4]

            lines!(ax3, points[:, 1], points[:, 2], points[:, 3], color = colors[i])

            if (cylindersInfo.cameraProjectionMatrix != undef)
                points2D = [cylindersInfo.cameraProjectionMatrix * point for point in eachrow(points)]
                points2D = [(point ./ point[3]) for point in points2D]
                points2D = hcat(points2D...)'

                lines!(ax2, points2D[:, 1], -points2D[:, 2], color = colors[i])
            end
        end
    end

    function plot2DPoints(singularPoints)
        for (i, singularPoint) in enumerate(singularPoints)
            scatter!(ax2, (singularPoint[1], -singularPoint[2]), color = colors[i])
        end
    end

    function plot2DCylinders(conicBorders; linestyle = :solid)
        y = function (x, l) return (-(l[1] * x + l[3]) / l[2]) end
        for (i, conicBorder) in enumerate(conicBorders)
            for line in conicBorder
                y1 = function (x) return y(x, line) end
                xs = -50:1:50
                ys1 = y1.(xs)
                lines!(ax2, xs, -ys1, color = colors[i], linestyle=linestyle)
            end
        end
    end
end