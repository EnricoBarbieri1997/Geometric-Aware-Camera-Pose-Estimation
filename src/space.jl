module Space
	export Vec3, Vec3Tuple, Point, PointTuple, transformation, identity_transformation, random_transformation, position_rotation, build_rotation_matrix

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

	function RotRad(x, y, z)
		return RotZYX(z, y, x)
	end

	function RotDeg(x, y, z)
		# Convert degrees to radians
		x = deg2rad(x)
		y = deg2rad(y)
		z = deg2rad(z)

		return RotRad(x, y, z)
	end

	function transformation(
		translation::Union{Array{<:Number, 1}, Vector{<:Number}} = [0, 0, 0],
		rotation::Union{Array{<:Number, 1}, Vector{<:Number}} = [0, 0, 0]
	)
		rotation = RotDeg(rotation...)

		r = zeros(4, 4)
		r[1:3, 1:3] .= (rotation)
		r[4, 4] = 1
		t = diagm([1.0, 1.0, 1.0, 1.0])
		t[1:3, 4] .= (translation .* 1)
	
		transform = t * r
		transform[abs.(transform) .< 1e-11] .= 0.0

		return transform
	end
	
	function position_rotation(transform_matrix::Matrix{<:Number})
		# Extract the translation vector from the transformation matrix
		translation = transform_matrix[1:3, 4]

		# Extract the rotation matrix from the transformation matrix
		rotation_matrix = transform_matrix[1:3, 1:3]

		# Convert the rotation matrix to Euler angles (in degrees)
		euler_angles = eulerangles_from_rotationmatrix(rotation_matrix)

		return translation, euler_angles
	end

	function random_transformation(
		centerBoundaries::Array{<:Number, 2} = [[-5, 5], [-5, 5], [-5, 5]]
	)
		center = rand_in_range(collect(centerBoundaries))
		rotation = rand_in_range((0, 40), 3)

		return transformation((center[1], center[2], center[3]), (rotation[1], rotation[2], rotation[3]))
	end

	function identity_transformation()
		return transformation()
	end

	function build_rotation_matrix(x, y, z, include_normalization = false)
		# R parametrized by x, y, z
		# https://en.wikipedia.org/wiki/Cayley_transform#Examples
		# 4.1.2 https://www.cv-foundation.org/openaccess/content_cvpr_2016/papers/Kukelova_Efficient_Intersection_of_CVPR_2016_paper.pdf
		k = 1 + x^2 + y^2 + z^2
		Rₚ = [
				1 + x^2 - y^2 - z^2     2*x*y - 2*z        2*y + 2*x*z;
				2*z + 2*x*y             1 - x^2 + y^2 - z^2  2*y*z - 2*x;
				2*x*z - 2*y             2*x + 2*y*z        1 - x^2 - y^2 + z^2
		]

		if (include_normalization)
				Rₚ = (1 / k) * Rₚ
		end

		return Rₚ
	end
end
