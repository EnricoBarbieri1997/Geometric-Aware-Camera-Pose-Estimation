module Cylinder
	export FromTransform, Random

	using ..Space, ..Utils, LinearAlgebra, Random, Rotations

	function FromTransform(center::PointTuple = (0, 0, 0), radius::Number = 1, rotation::Vec3Tuple = (0, 0, 0))
		rotation = RotXYZ(deg2rad.(rotation)...)
		canonicalCylinder = diagm([1, 1, 0, -(radius^2)])

		r = zeros(4, 4)
		r[1:3, 1:3] .= (rotation)
		r[4, 4] = 1
		t = diagm([1.0, 1.0, 1.0, 1.0])
		t[1:3, 4] .= (center .* -1)
		
		transformMatrix = r * t

		display(transformMatrix)
		
		cylinder = transformMatrix' * canonicalCylinder * transformMatrix
		return cylinder
	end

	function Random(
		centerBoundaries::Tuple{Tuple{Number, Number}, Tuple{Number, Number}, Tuple{Number, Number}} = ((-5, 5), (-5, 5), (-5, 5)),
		radiusBoundaries::Tuple{Number, Number} = (1, 3),
	)
		center = (randRange(centerBoundaries[1]), randRange(centerBoundaries[2]), randRange(centerBoundaries[3]))
		radius = randRange(radiusBoundaries)
		rotation = (randRange((0, 360)), randRange((0, 360)), randRange((0, 360)))
		return FromTransform(center, radius, rotation)
	end
end