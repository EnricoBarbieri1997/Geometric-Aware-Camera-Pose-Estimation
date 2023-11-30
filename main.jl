include("includes.jl")

using .Space: Transformation, RandomTransformation
using .Utils
using LinearAlgebra: deg2rad, diagm, normalize, svd
using GLMakie, GLMakie.FileIO, HomotopyContinuation, Polynomials, Rotations

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

    function lines_from_conic(i)
        @var x y z
        quadric = cylinders[i][2]
        f₁ = quadric[1] * x^2 + 2 * quadric[2] * x * y + 2 * quadric[3] * x * z + 2 * quadric[4] * x + quadric[6] * y^2 + 2 * quadric[7] * y * z + 2 * quadric[8] * y + quadric[11] * z^2 + 2 * quadric[12] * z + quadric[16]
        f₂ = x * cameraOrigin[1] + y * cameraOrigin[2] + z * cameraOrigin[3] + cameraOrigin[4]
        f₃ = z - 1
        F = System([f₁, f₂, f₃])
        result = solve(F)
        planes = real_solutions(result)
        lines = [cameraProjectionMatrix * [plane... 1]' for plane in planes]
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