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
end

# cameraTranslation = (2.0, 30.0, 5.0)
# cameraRotation = (-83.0, 180.0, 0.0)
cameraTranslation = (0.0, 0.0, 0.0)
cameraRotation = (0.0, 0.0, 0.0)
cameraPositionMatrix = Transformation(cameraTranslation, cameraRotation)
cameraProjectionMatrix = CameraMatrix(cameraTranslation, cameraRotation, 1, 1)

conics = Array{Tuple{Matrix{Float64}, Tuple{Matrix{Float64}, Vector{Float64}}}}(undef, numberOfCylinders)
for i in 1:numberOfCylinders
	conics[i] = (cameraProjectionMatrix * cylinders[i][1] * cameraProjectionMatrix', (cameraProjectionMatrix * cylinders[i][2][1] * cameraProjectionMatrix', cameraProjectionMatrix * cylinders[i][2][2]))
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

display("Camera projection matrix: $(cameraProjectionMatrix)")
println("Singular Points:")
for i in 1:numberOfCylinders
    println("Original Cylinder $i Vanishing Point: ", cylinders[i][2][2])
    println("Cylinder $i: ", singularPoints[i])
end

# f = initFigure()
# plot3DCamera(Plot3DCameraInput(
#     cameraRotation,
#     cameraTranslation
# ))
# plot3DCylinders(Plot3DCylindersInput(
#     transforms,
#     radiuses,
#     numberOfCylinders,
#     # cameraProjectionMatrix
# ))
# plot2DPoints(singularPoints)
# plot2DCylinders(conicBorders)
# f

# 3 line minimum to solve the pose
# numberOfLinesToSolveFor = 3
# lines = []
# pointAtInfinityToUse = []
# dualQuadricToUse = []
# for (i, borders) in enumerate(conicBorders)
#     for line in borders
#         push!(lines, line)
#         push!(pointAtInfinityToUse, cylinders[i][2][2][1:3])
#         push!(dualQuadricToUse, cylinders[i][2][1])
#         if (size(lines)[1] == numberOfLinesToSolveFor)
#             break
#         end
#     end
#     if (size(lines)[1] == numberOfLinesToSolveFor)
#         break
#     end
# end

# return
# @var x y z
# # R parametrized by x, y, z
# # https://en.wikipedia.org/wiki/Cayley_transform#Examples
# # 4.1.2 https://www.cv-foundation.org/openaccess/content_cvpr_2016/papers/Kukelova_Efficient_Intersection_of_CVPR_2016_paper.pdf
# k = 1 + x^2 + y^2 + z^2
# Rₚ = #= (1/k) * =# [
#     1 + x^2 - y^2 - z^2     2*x*y - 2*z        2*y + 2*x*z;
#     2*z + 2*x*y             1 - x^2 + y^2 - z^2  2*y*z - 2*x;
#     2*x*z - 2*y             2*x + 2*y*z        1 - x^2 - y^2 + z^2
# ]

# systemToSolve = []
# for i in 1:numberOfLinesToSolveFor
#     equation = lines[i]' * Rₚ * pointAtInfinityToUse[i]
#     push!(systemToSolve, equation)
# end

# F = System(systemToSolve)
# result = solve(F)
# display(result)
# solution = real_solutions(result)[1]
# xₛ = solution[1]
# yₛ = solution[2]
# zₛ = solution[3]
# # Rotation as quaternion
# rotationCalculated = Rotations.QuatRotation(1, xₛ , yₛ, zₛ)
# display("Rotation angle: $(rotation_angle(rotationCalculated)), Rotation axis: $(rotation_axis(rotationCalculated))")

# @var tx ty tz
# P = [1 0 0 0;
#     0 1 0 0;
#     0 0 1 0] * [rotationCalculated [tx; ty; tz]; 0 0 0 1]

# systemToSolve = []
# for i in 1:3
#     equation = lines[i]' * P * dualQuadricToUse[i] * P' * lines[i]
#     push!(systemToSolve, equation)
# end

# F = System(systemToSolve)
# result = solve(F)
# solution = real_solutions(result)[1]
# translationCalculated = solution
# display("Translation: $(translationCalculated)")

# calculatedCameraMatrix = CameraMatrix((translationCalculated[1], translationCalculated[2], translationCalculated[3]), (xₛ, yₛ, zₛ), 1, 1)

# display("Camera projection matrix: $(cameraProjectionMatrix ./ cameraProjectionMatrix[3, 4])")
# display("Calculated projection camera matrix: $(calculatedCameraMatrix ./ calculatedCameraMatrix[3, 4])")