include("includes.jl")

using .Space: Transformation

numberOfCylinders = 4
cylinders = Array{Tuple{Matrix{Float64}, Matrix{Float64}}}(undef, numberOfCylinders)
for i in 1:numberOfCylinders
	cylinders[i] = Cylinder.StandardAndDualRandom()

	display(Quadric.ToFormula(cylinders[i][1]))
end

cameraMatrix = Transformation((1, 3, 7), (-40, 0, -15))
cameraMatrix = cameraMatrix ./ cameraMatrix[4, 4]
cameraProjectionMatrix = cameraMatrix[1:3, :]

conics = Array{Matrix{Float64}}(undef, numberOfCylinders)
dualConics = Array{Matrix{Float64}}(undef, numberOfCylinders)
for i in 1:numberOfCylinders
	conics[i] = cameraProjectionMatrix * cylinders[i][1] * cameraProjectionMatrix'
	dualConics[i] = cameraProjectionMatrix * cylinders[i][2] * cameraProjectionMatrix'

	# display(Conic.ToFormula(conics[i]))
end

