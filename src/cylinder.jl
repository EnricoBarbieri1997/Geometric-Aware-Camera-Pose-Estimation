module Cylinder
	export CylinderProperties, standard_and_dual, standard, dual, standard_and_dual_random, standard_random, random_dual

	using ..Space, ..Utils, LinearAlgebra, Random

	mutable struct CylinderProperties
		euler_rotation::Vector{Number}
		radiuses::Array{<:Number, 1}
		matrix::Matrix{<:Number}
		singular_point::Vector{Number}
		dual_matrix::Matrix{<:Number}
		transform::Matrix{<:Number}

		CylinderProperties() = new()
	end

	function standard_and_dual(
		transform_matrix::Matrix{<:Number},
		radius::Union{Vector{<:Number}, Array{<:Number, 1}} = [1, 1]
	)
		inverseRadiusSquareX = 1 / (radius[1]^2)
		inverseRadiusSquareY = 1 / (radius[2]^2)
		canonicalCylinder = diagm([inverseRadiusSquareX, inverseRadiusSquareY, 0, -1])

		cylinder = transform_matrix' * canonicalCylinder * transform_matrix

		dualCanonicalCylinderMatrix = zeros(4, 4)
		dualCanonicalCylinderMatrix[[1, 2, 4], [1, 2, 4]] .= inv(canonicalCylinder[[1, 2, 4], [1, 2, 4]])

		dualCylinderMatrix = inv(transform_matrix) * dualCanonicalCylinderMatrix * inv(transform_matrix')
		dualCylinderSingularPoint = transform_matrix * [0; 0; 1; 0]

		return cylinder, dualCylinderMatrix, dualCylinderSingularPoint
	end

	function standard_and_dual(
		center::Vector{<:Number} = [0, 0, 0],
		radius::Vector{<:Number} = [1, 1],
		rotation::Vector{<:Number} = [0, 0, 0]
	)
		transform_matrix = transformation(center, rotation)
		return standard_and_dual(transform_matrix, radius)
	end

	function standard_and_dual(
		center::Vector{<:Number} = [0, 0, 0],
		radius::Number = 1,
		rotation::Vector{<:Number} = [0, 0, 0]
	)
		return standard_and_dual(center, (radius, radius), rotation)
	end

	function standard(
		center::Vector{<:Number} = [0, 0, 0],
		radius::Number = 1,
		rotation::Vector{<:Number} = [0, 0, 0]
	)
		return standard_and_dual(center, radius, rotation)[1]
	end

	function dual(
		center::Vector{<:Number} = [0, 0, 0],
		radius::Number = 1,
		rotation::Vector{<:Number} = [0, 0, 0]
	)
		return standard_and_dual(center, radius, rotation)[2]
	end

	function standard_and_dual_random(
		centerBoundaries::Tuple{Tuple{Number, Number}, Tuple{Number, Number}, Tuple{Number, Number}} = ((-5, 5), (-5, 5), (-5, 5)),
		radiusBoundaries::Tuple{Number, Number} = (1, 3),
	)
		center = rand_in_range(collect(centerBoundaries))
		radius = rand_in_range(radiusBoundaries, 2)
		rotation = rand_in_range((0, 360), 3)
		return standard_and_dual((center[1], center[2], center[3]), (radius[1], radius[2]), (rotation[1], rotation[2], rotation[3]))
	end

	function standard_random(
		centerBoundaries::Tuple{Tuple{Number, Number}, Tuple{Number, Number}, Tuple{Number, Number}} = ((-5, 5), (-5, 5), (-5, 5)),
		radiusBoundaries::Tuple{Number, Number} = (1, 3),
	)
		return standard_and_dual_random(centerBoundaries, radiusBoundaries)[1]
	end

	function random_dual(
		centerBoundaries::Tuple{Tuple{Number, Number}, Tuple{Number, Number}, Tuple{Number, Number}} = ((-5, 5), (-5, 5), (-5, 5)),
		radiusBoundaries::Tuple{Number, Number} = (1, 3),
	)
		return standard_and_dual_random(centerBoundaries, radiusBoundaries)[2]
	end
end