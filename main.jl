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

cameraTranslation = (2.0, 30.0, 5.0)
cameraRotation = (-83.0, 180.0, 0.0)
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
lines = []
pointAtInfinityToUse = []
dualQuadricToUse = []
for (i, borders) in enumerate(conicBorders)
    for line in borders
        push!(lines, line)
        push!(pointAtInfinityToUse, cylinders[i][2][2][1:3])
        push!(dualQuadricToUse, cylinders[i][2][1])
        if (size(lines)[1] == 3)
            break
        end
    end
    if (size(lines)[1] == 3)
        break
    end
end

return
@var α β γ
Rx = [1 0 0; 0 cos(α) -sin(α); 0 sin(α) cos(α)]
Ry = [cos(β) 0 sin(β); 0 1 0; -sin(β) 0 cos(β)]
Rz = [cos(γ) -sin(γ) 0; sin(γ) cos(γ) 0; 0 0 1]

R = Rx * Ry * Rz

systemToSolve = []
for i in 1:3
    equation = lines[i]' * R * pointAtInfinityToUse[i]
    display(equation)
    display("\n")
    push!(systemToSolve, equation)
end

F = System(systemToSolve)
result = solve(F)
# display(result)
# solution = real_solutions(result)[1]
# αₛ = solution[1]
# βₛ = solution[2]
# γₛ = solution[3]
# Rx = [1 0 0; 0 cos(αₛ) -sin(αₛ); 0 sin(αₛ) cos(αₛ)]
# Ry = [cos(βₛ) 0 sin(βₛ); 0 1 0; -sin(βₛ) 0 cos(βₛ)]
# Rz = [cos(γₛ) -sin(γₛ) 0; sin(γₛ) cos(γₛ) 0; 0 0 1]
# R = Rx * Ry * Rz

# @var tx ty tz
# P = [R [tx; ty; tz]; 0 0 0 1]

# systemToSolve = []
# for i in 1:3
#     equation = lines[i]' * P * dualQuadricToUse[i] * P' * lines[i]
#     push!(systemToSolve, equation)
# end

# F = System(systemToSolve)
# result = solve(F)