include("includes.jl")

using .Space: Transformation, RandomTransformation, IdentityTransformation
using .Utils
using LinearAlgebra: deg2rad, diagm, dot, normalize, svd
using GLMakie, GLMakie.FileIO, HomotopyContinuation, Polynomials, Rotations

numberOfCylinders = 4
cylinders = Array{Tuple{Matrix{Float64}, Tuple{Matrix{Float64}, Vector{Float64}}}}(undef, numberOfCylinders)
transforms = Array{Matrix{Float64}}(undef, numberOfCylinders)
radiuses = Array{Tuple{Number, Number}}(undef, numberOfCylinders)
for i in 1:numberOfCylinders
    transforms[i] = RandomTransformation()
    radius = randRange((1, 3), 2)
    radiuses[i] = (radius[1], radius[2])
	cylinders[i] = Cylinder.StandardAndDual(transforms[i], radiuses[i])

	display(Quadric.ToFormula(cylinders[i][1]))
    display(PointFormulas.ToFormula(cylinders[i][2][2] ./ cylinders[i][2][2][4]))
end

cameraTranslation = (2.0, 30.0, 5.0)
cameraRotation = (-83.0, 0.0, 180.0)
cameraMatrix = Transformation(cameraTranslation, cameraRotation)
cameraMatrix = cameraMatrix ./ cameraMatrix[4, 4]
# cameraMatrix = Transformation((0,0,0), (0, 90, 0)) * cameraMatrix
focalLength = 1.0
pinHolePerfectModel = [1 0 0 0; 0 1 0 0; 0 0 1/focalLength 0]
cameraOrigin = cameraMatrix * [0.0, 0.0, 0.0, 1.0]
cameraOrigin = cameraOrigin ./ cameraOrigin[4]
# display(PointFormulas.ToFormula(cameraOrigin))
cameraProjectionMatrix = pinHolePerfectModel * cameraMatrix
cameraProjectionMatrix = cameraProjectionMatrix ./ cameraProjectionMatrix[3, 4]
relativeZero = 100
cameraProjectionMatrix = floor.(cameraProjectionMatrix .* relativeZero) ./ relativeZero

conics = Array{Tuple{Matrix{Float64}, Tuple{Matrix{Float64}, Vector{Float64}}}}(undef, numberOfCylinders)
for i in 1:numberOfCylinders
	conics[i] = (cameraProjectionMatrix * cylinders[i][1] * cameraProjectionMatrix', (cameraProjectionMatrix * cylinders[i][2][1] * cameraProjectionMatrix', cameraProjectionMatrix * cylinders[i][2][2]))
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

    z, θ = LinRange(-20, 20, heightLevels), LinRange(0, 2π, angles)
    x = cos.(θ)
    y = sin.(θ)

    for i in 1:numberOfCylinders
        radius = radiuses[i]
        X = radius[1] * x
        Y = radius[2] * y

        canonicPoints = []

        for j in 1:heightLevels
            canonicPoints = vcat(canonicPoints, [X Y (z[j] * ones(angles)) ones(angles)])
        end
        points = transpose(transforms[i] * canonicPoints')
        points = points ./ points[:, 4]
        
        lines!(ax3, points[:, 1], points[:, 2], points[:, 3], color = colors[i])

        points2D = [cameraProjectionMatrix * point for point in eachrow(points)]
        # points2D = [(point ./ point[3]) for point in points2D]
        points2D = hcat(points2D...)'

        lines!(ax2, points2D[:, 1], points2D[:, 2], color = colors[i])
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
        line = [x, y, z]
        conicSingularPoint = conics[i][2][2]
        conicMatrix = conics[i][2][1]
        f₁ = line' * conicSingularPoint
        f₂ = line' * conicMatrix * line
        f₃ = z - 1
        F = System([f₁, f₂, f₃])
        result = solve(F)
        lines = result
        return real_solutions(lines)
    end

    for i in 1:numberOfCylinders
        singularPoint = conics[i][2][2]
        singularPoint = singularPoint ./ singularPoint[3]
        # singularPoint = singularPoint .* 500
        singularPoint = (singularPoint[1], singularPoint[2])
        scatter!(singularPoint, color = colors[i])
        # lines = lines_from_conic(i)
        # y = function (x, l) return (-(l[1] * x + l[3]) / l[2]) end
        # for line in lines
        #     y1 = function (x) return y(x, line) end
        #     xs = -10:1:10
        #     ys1 = y1.(xs)
        #     lines!(ax2, xs, ys1, color = colors[i])
        # end
    end
end
plot2D()

f