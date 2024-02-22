module Space
	export Vec3, Vec3Tuple, Point, PointTuple, Transformation, IdentityTransformation, RandomTransformation

	using ..Utils
	using LinearAlgebra
	using Rotations

	struct Vec3
		x::Number
		y::Number
		z::Number
	end
	const Vec3Tuple = Tuple{Number, Number, Number}

	const Point = Vec3
	const PointTuple = Vec3Tuple

	function Transformation(
		translation::Vec3Tuple = (0, 0, 0),
		rotation::Vec3Tuple = (0, 0, 0)
	)
		rotation = RotXYZ(deg2rad.(rotation)...)

		r = zeros(4, 4)
		r[1:3, 1:3] .= (rotation)
		r[4, 4] = 1
		t = diagm([1.0, 1.0, 1.0, 1.0])
		t[1:3, 4] .= (translation .* 1)

		return t * r
	end

	function RandomTransformation(
		centerBoundaries::Tuple{Tuple{Number, Number}, Tuple{Number, Number}, Tuple{Number, Number}} = ((-5, 5), (-5, 5), (-5, 5))
	)
		center = randRange(collect(centerBoundaries))
		rotation = randRange((0, 360), 3)

		return Transformation((center[1], center[2], center[3]), (rotation[1], rotation[2], rotation[3]))
	end

	function IdentityTransformation()
		return Transformation()
	end
end
