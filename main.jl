include("includes.jl")

using .Space: Transformation, RandomTransformation, IdentityTransformation
using .Camera: CameraMatrix
using .Plotting: initFigure, plot2DPoints, Plot3DCameraInput, plot3DCamera, Plot3DCylindersInput, plot3DCylinders, plot2DCylinders
using .Debug
using .Utils
using LinearAlgebra: deg2rad, diagm, dot, normalize, svd
using HomotopyContinuation, Polynomials, Rotations

numberOfCylinders = 4
cylinders = Array{Tuple{Matrix{Float64}, Tuple{Matrix{Float64}, Vector{Float64}}}}(undef, numberOfCylinders)
transforms = Array{Matrix{Float64}}(undef, numberOfCylinders)
radiuses = Array{Tuple{Number, Number}}(undef, numberOfCylinders)
for i in 1:numberOfCylinders
    transforms[i] = RandomTransformation()
    radius = randRange((1, 3), 2)
    radiuses[i] = (radius[1], radius[2])
	cylinders[i] = Cylinder.StandardAndDual(transforms[i], radiuses[i])

	# display(Quadric.ToFormula(cylinders[i][1]))
    # display(PointFormulas.ToFormula(cylinders[i][2][2] ./ cylinders[i][2][2][4]))
end

cameraTranslation = (2.0, 30.0, 5.0)
cameraRotation = (-83.0, 180.0, 0.0)
cameraPositionMatrix = Transformation(cameraTranslation, cameraRotation)
cameraPositionMatrix = cameraPositionMatrix ./ cameraPositionMatrix[4, 4]
# cameraPositionMatrix = Transformation((0,0,0), (0, 90, 0)) * cameraPositionMatrix
cameraProjectionMatrix = CameraMatrix(cameraTranslation, cameraRotation, 72, 0.55)
cameraOrigin = cameraPositionMatrix * [0.0, 0.0, 0.0, 1.0]
cameraOrigin = cameraOrigin ./ cameraOrigin[4]
# display(PointFormulas.ToFormula(cameraOrigin))
relativeZero = 100
cameraProjectionMatrix = floor.(cameraProjectionMatrix .* relativeZero) ./ relativeZero

conics = Array{Tuple{Matrix{Float64}, Tuple{Matrix{Float64}, Vector{Float64}}}}(undef, numberOfCylinders)
for i in 1:numberOfCylinders
	conics[i] = (cameraProjectionMatrix * cylinders[i][1] * cameraProjectionMatrix', (cameraProjectionMatrix * cylinders[i][2][1] * cameraProjectionMatrix', cameraProjectionMatrix * cylinders[i][2][2]))
	# display(Conic.ToFormula(conics[i][1]))
end

singularPoints = Array{Tuple{Number, Number}}(undef, numberOfCylinders)

for i in 1:numberOfCylinders
    singularPoint = conics[i][2][2]
    singularPoint = singularPoint ./ singularPoint[3]
    # singularPoint = singularPoint .* 500
    singularPoint = (singularPoint[1], singularPoint[2])
    singularPoints[i] = singularPoint
end

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

conicBorders = Array{Array{Vector{Float64}}}(undef, numberOfCylinders)
for i in 1:numberOfCylinders
    lines = lines_from_conic(i)
    conicBorders[i] = Array{Vector{Float64}}(undef, length(lines))
    for (j, line) in enumerate(lines)
        conicBorders[i][j] = line
    end
end

f = initFigure()
plot3DCamera(Plot3DCameraInput(
    cameraRotation,
    cameraTranslation
))
plot3DCylinders(Plot3DCylindersInput(
    transforms,
    radiuses,
    numberOfCylinders,
    # cameraProjectionMatrix
))
plot2DPoints(singularPoints)
plot2DCylinders(conicBorders)
f