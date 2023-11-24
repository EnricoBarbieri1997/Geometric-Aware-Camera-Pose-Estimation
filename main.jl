include("includes.jl")

numberOfCylinders = 4
cylinders = Array{Matrix{Float64}}(undef, numberOfCylinders)
for i in 1:numberOfCylinders
	cylinders[i] = Cylinder.Random()
end

Quadric.ToFormula(cylinders[1])