include("includes.jl")

using .Space: Transformation, RandomTransformation
using .Utils
using LinearAlgebra: deg2rad, diagm, normalize, svd
using GLMakie, GLMakie.FileIO, Polynomials, Rotations

numberOfCylinders = 4
cylinders = Array{Tuple{Matrix{Float64}, Matrix{Float64}}}(undef, numberOfCylinders)
transforms = Array{Matrix{Float64}}(undef, numberOfCylinders)
radiuses = Array{Tuple{Number, Number}}(undef, numberOfCylinders)
for i in 1:numberOfCylinders
    transforms[i] = RandomTransformation()
    radius = randRange((1, 3), 2)
    radiuses[i] = (radius[1], radius[2])
	cylinders[i] = Cylinder.StandardAndDual(transforms[i], radiuses[i])

	# display(Quadric.ToFormula(cylinders[i][1]))
end

cameraTranslation = (2.0, 30.0, 5.0)
cameraRotation = (-83.0, 0.0, 180.0)
cameraMatrix = Transformation(cameraTranslation, cameraRotation)
cameraMatrix = cameraMatrix ./ cameraMatrix[4, 4]
cameraOrigin = cameraMatrix * [0.0, 0.0, 0.0, 1.0]
cameraOrigin = cameraOrigin ./ cameraOrigin[4]
# display(PointFormulas.ToFormula(cameraOrigin))
cameraProjectionMatrix = cameraMatrix[1:3, :]

conics = Array{Tuple{Matrix{Float64}, Matrix{Float64}}}(undef, numberOfCylinders)
for i in 1:numberOfCylinders
	conics[i] = (cameraProjectionMatrix * cylinders[i][1] * cameraProjectionMatrix', cameraProjectionMatrix * cylinders[i][2] * cameraProjectionMatrix')

	# display(Conic.ToFormula(conics[i][1]))
end

f = Figure(size=(1200, 800))
ax3 = Axis3(f[1, 1], title = "Cylinders", aspect = :equal)
ax2 = Axis(f[1, 2], title = "Conics", autolimitaspect = 1)

colors = [:red, :green, :blue, :yellow]

points3d = undef

function plot3D()
    global colors, cameraRotation, cameraTranslation, transforms, radiuses, numberOfCylinders, ax3
    heightLevels = 100
    angles = 100

    z, θ = LinRange(-10, 20, heightLevels), LinRange(0, 2π, angles)
    x = cos.(θ)
    y = sin.(θ)

    for i in 1:numberOfCylinders
        radius = radiuses[i]
        X = radius[1] * x
        Y = radius[2] * y

        for j in 1:heightLevels
            points = [X Y (z[j] * ones(angles)) ones(angles)] * transforms[i]'
            points = points ./ points[:, 4]
            
            lines!(ax3, points[:, 1], points[:, 2], points[:, 3], color = colors[i])

            points2D = [cameraProjectionMatrix * point for point in eachrow(points)]
            points2D = vcat(transpose.(points2D)...)

            lines!(ax2, points2D[:, 1], points2D[:, 2], color = colors[i])
        end
    end

    camera = load("./camera.stl")
    cameraMesh = mesh!(
        ax3,
        camera,
    )
    cameraRotationRad = deg2rad.(cameraRotation)
    cameraRotation = RotXYZ(cameraRotationRad...)
    cameraRotationAxis = rotation_axis(cameraRotation)
    cameraRotationAngle = rotation_angle(cameraRotation)
    rotate!(cameraMesh, cameraRotationAxis, cameraRotationAngle)
    translate!(cameraMesh, cameraTranslation)
end
plot3D()

function plot2D()
    global colors, radiuses, transforms, cameraProjectionMatrix, numberOfCylinders, ax2
    θ = LinRange(0, 2π, 100)

    function lines_from_conic(i)
        planes = []
        # display(i)
        for angle in θ
            x = radiuses[i][1] * cos(angle)
            y = radiuses[i][2] * sin(angle)
            planeDirection = normalize([x, y, 0])
            plane = [planeDirection... (-sqrt(x^2 + y^2))]
            plane = transforms[i]' * plane'
            distance = (cameraOrigin' * plane)[1]
            if(abs(distance) <= 1)
                # display(Plane.ToFormula(plane))
                push!(planes, plane)
            end
        end
        lines = [cameraProjectionMatrix * plane for plane in planes]
        return lines
    end

    for i in 1:numberOfCylinders
        lines = lines_from_conic(i)
        y = function (x, l) return (-(l[1] * x + l[3]) / l[2]) end
        for line in lines
            y1 = function (x) return y(x, line) end
            xs = -10:1:10
            ys1 = y1.(xs)
            lines!(ax2, xs, ys1, color = colors[i])
        end
    end
end
# plot2D()

f