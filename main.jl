include("includes.jl")

using .Space: Transformation, RandomTransformation
using .Utils
using LinearAlgebra: deg2rad, diagm, svd
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
cameraRotation = (7.0, 0.0, 180.0)
cameraMatrix = Transformation(cameraTranslation, cameraRotation)
cameraMatrix = cameraMatrix ./ cameraMatrix[4, 4]
cameraProjectionMatrix = cameraMatrix[1:3, :]

conics = Array{Tuple{Matrix{Float64}, Matrix{Float64}}}(undef, numberOfCylinders)
for i in 1:numberOfCylinders
	conics[i] = (cameraProjectionMatrix * cylinders[i][1] * cameraProjectionMatrix', cameraProjectionMatrix * cylinders[i][2] * cameraProjectionMatrix')

	# display(Conic.ToFormula(conics[i][1]))
end

f = Figure(size=(1200, 800))
ax3 = Axis3(f[1, 1], title = "Cylinders", aspect = (1, 3, 1))
ax2 = Axis(f[1, 2], title = "Conics", aspect = 1)

colors = [:red, :green, :blue, :yellow]

function plot3D()
    global colors, cameraRotation, cameraTranslation, transforms, radiuses, numberOfCylinders, ax3
    heightLevels = 100
    angles = 100

    z, θ = LinRange(0, 10, heightLevels), LinRange(0, 2π, angles)
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
    global colors
    function lines_from_conic(dualConic, originalTransformation, cameraMatrix)
        projectedZPlane = cameraMatrix * (inv(originalTransformation) * diagm([0, 0, 1, 0]) * inv(originalTransformation')) * cameraMatrix'
        a, b, v = svd(projectedZPlane)
        b1 = v[:, 2]
        b2 = v[:, 3]
        cc = [b1'*dualConic*b1 2*b1'*dualConic*b2 b2'*dualConic*b2]
        ll = roots(Polynomial(cc))
        l1 = b1*ll[1]+b2
        l2 = b1*ll[2]+b2
        l1 = l1 ./ l1[3]
        l2 = l2 ./ l2[3]
        return l1, l2
    end

    for i in 1:numberOfCylinders
        l1, l2 = lines_from_conic(conics[i][2], transforms[i], cameraProjectionMatrix)
        y = function (x, l) return (-(l[1] * x + l[3]) / l[2]) end
        y1 = function (x) return y(x, l1) end
        y2 = function (x) return y(x, l2) end
        xs = -10:20:10
        ys1 = y1.(xs)
        ys2 = y2.(xs)
        lines!(ax2, xs, ys1, color = colors[i])
        lines!(ax2, xs, ys2, color = colors[i])
    end
end
plot2D()

f