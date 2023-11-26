module Cylinder
	export StandardAndDual, Standard, Dual, StandardAndDualRandom, StandardRandom, RandomDual

	using ..Space, ..Utils, LinearAlgebra, Random

	function StandardAndDual(
		transformMatrix::Matrix{Float64},
		radius::Tuple{Number, Number} = (1, 1)
	)
		inverseRadiusSquareX = 1 / (radius[1]^2)
		inverseRadiusSquareY = 1 / (radius[2]^2)
		canonicalCylinder = diagm([inverseRadiusSquareX, inverseRadiusSquareY, 0, -1])

		cylinder = transformMatrix' * canonicalCylinder * transformMatrix

		dualCanonicalCylinder = zeros(4, 4)
		dualCanonicalCylinder[[1, 2, 4], [1, 2, 4]] .= inv(canonicalCylinder[[1, 2, 4], [1, 2, 4]])

		dualCylinder = inv(transformMatrix) * dualCanonicalCylinder * inv(transformMatrix')
		return cylinder, dualCylinder
	end

	function StandardAndDual(
		center::PointTuple = (0, 0, 0),
		radius::Tuple{Number, Number} = (1, 1), rotation::Vec3Tuple = (0, 0, 0)
	)
		transformMatrix = Transformation(center, rotation)
		return StandardAndDual(transformMatrix, radius)
	end

	function StandardAndDual(
		center::PointTuple = (0, 0, 0),
		radius::Number = 1, rotation::Vec3Tuple = (0, 0, 0)
	)
		return StandardAndDual(center, (radius, radius), rotation)
	end

	function Standard(
		center::PointTuple = (0, 0, 0),
		radius::Number = 1, rotation::Vec3Tuple = (0, 0, 0)
	)
		return StandardAndDual(center, radius, rotation)[1]
	end

	function Dual(
		center::PointTuple = (0, 0, 0),
		radius::Number = 1, rotation::Vec3Tuple = (0, 0, 0)
	)
		return StandardAndDual(center, radius, rotation)[2]
	end

	function StandardAndDualRandom(
		centerBoundaries::Tuple{Tuple{Number, Number}, Tuple{Number, Number}, Tuple{Number, Number}} = ((-5, 5), (-5, 5), (-5, 5)),
		radiusBoundaries::Tuple{Number, Number} = (1, 3),
	)
		center = randRange(collect(centerBoundaries))
		radius = randRange(radiusBoundaries, 2)
		rotation = randRange((0, 360), 3)
		return StandardAndDual((center[1], center[2], center[3]), (radius[1], radius[2]), (rotation[1], rotation[2], rotation[3]))
	end

	function StandardRandom(
		centerBoundaries::Tuple{Tuple{Number, Number}, Tuple{Number, Number}, Tuple{Number, Number}} = ((-5, 5), (-5, 5), (-5, 5)),
		radiusBoundaries::Tuple{Number, Number} = (1, 3),
	)
		return StandardAndDualRandom(centerBoundaries, radiusBoundaries)[1]
	end

	function RandomDual(
		centerBoundaries::Tuple{Tuple{Number, Number}, Tuple{Number, Number}, Tuple{Number, Number}} = ((-5, 5), (-5, 5), (-5, 5)),
		radiusBoundaries::Tuple{Number, Number} = (1, 3),
	)
		return StandardAndDualRandom(centerBoundaries, radiusBoundaries)[2]
	end
end