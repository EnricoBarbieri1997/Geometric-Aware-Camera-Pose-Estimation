module Geometry
	export Plane, Line, Point, Circle, TangentLineNotFound, homogeneous_line_from_points, issame_line, rotation_between_lines, cylinder_rotation_from_axis, homogeneous_to_line, line_to_homogenous, homogeneous_line_intercept, homogeneous_anglebetween, project_point_into_line, project_point_into_plane, get_tangentpoints_circle_point, get_cylinder_contours

	using ..Utils
	using ..Camera: CameraProperties
	using ..Cylinder: CylinderProperties
	using LinearAlgebra: Diagonal, cross, dot, norm, normalize, pinv, svd
	using Rotations: RotMatrix, AngleAxis, RotXYZ, params as rotations_params
	using Polynomials: Polynomial, roots

	struct Point
		x::Number
		y::Number
	end
	struct Line
		origin::Vector{<:Number}
		direction::Vector{<:Number}
	end
	struct Plane
		origin::Vector{<:Number}
		normal::Vector{<:Number}
	end
	struct Circle
		center::Vector{<:Number}
		radius::Number
		axis::Union{Vector{<:Number}, Nothing}
	end

	struct TangentLineNotFound <: Exception
		msg::String
	end

	function homogeneous_line_from_points(p1, p2)
			# Convert to homogeneous coordinates
			p1h = [p1[1], p1[2], 1.0]
			p2h = [p2[1], p2[2], 1.0]

			# Line is the cross product of the two points
			l = cross(p1h, p2h)

			return l  # line coefficients [a, b, c] such that ax + by + c = 0
	end

	function issame_line(line₁::Line, line₂::Line, digits=6)
		direction_factors = line₁.direction ./ line₂.direction
		direction_factors = normalize(direction_factors)
		direction_factors = round.(direction_factors; digits=digits)
		issame_direction = allequal(direction_factors)
		if !issame_direction
			return false
		end

		origin_factors = (line₁.origin - line₂.origin) ./ line₁.direction
		origin_factors = normalize(origin_factors)
		origin_factors = round.(origin_factors; digits=digits)
		return allequal(origin_factors)
	end

	# x*l[1] + y*l[2] + l[3] = 0
	# y = -l[1]/l[2] * x - l[3]/l[2]
	function homogeneous_to_line(homogenous_line::Vector{<:Number})
		m = - homogenous_line[1] / homogenous_line[2]
		q = - homogenous_line[3] / homogenous_line[2]
		origin = [0, q]
		direction = [1, m + q] - origin
		return Line(origin, direction)
	end

	function line_to_homogenous(line::Line)
		homogenous_line = cross([line.origin; 1], [line.origin; 1] + [line.direction*10; 0])
		return homogenous_line ./ homogenous_line[3]
	end

	function homogeneous_line_intercept(x, line)
		return -line[1]/line[2] * x - line[3]/line[2]
	end

	function homogeneous_anglebetween(a, b)
		return atan((b[1]*a[2]-a[1]*b[2])/(a[1]*b[1]+a[2]*b[2]))
	end

	function project_point_into_line(point::Vector{<:Number}, line::Line)::Vector{<:Number}
		direction = line.direction / norm(line.direction)
		origin = line.origin
		v = point - origin
		return origin + dot(v, direction) * direction
	end

	function project_point_into_plane(point::Vector{<:Number}, plane::Plane)
		normal = plane.normal
		origin = plane.origin
		v = point - origin
		d = dot(v, normal)
		return point - d * normal
	end

	function rotation_between_lines(line₁::Vector{Float64}, line₂::Vector{Float64})
		# vertical axis (source)
		v1 = normalize(line₁)

		# your target axis (example)
		v2 = normalize(line₂)  # make sure it's normalized

		# compute rotation axis: cross product
		axis = cross(v1, v2)

		# compute angle: arccos of dot product
		angle = acos(clamp(dot(v1, v2), -1.0, 1.0))

		# handle the case when vectors are parallel or anti-parallel
		if isapprox(norm(axis), 0.0)
				if dot(v1, v2) > 0
						rot = one(RotMatrix{3})
				else
						# 180-degree rotation around any perpendicular axis
						rot = AngleAxis(π, 1.0, 0.0, 0.0)  # choose X axis arbitrarily
				end
		else
				axis = normalize(axis)
				rot = AngleAxis(angle, axis...)
		end

		# now rot is the rotation you want

		# if you want Euler angles (for example RotXYZ convention)
		euler = RotXYZ(rot)
		return rotations_params(euler)
	end

	function plane_through_3_points(p1::Vector{<:Number}, p2::Vector{<:Number}, p3::Vector{<:Number})
		# Create a plane through three points
		v1 = p2 - p1
		v2 = p3 - p1
		normal = cross(v1, v2)
		return Plane(p1, normal)
	end

	function plane_to_homogeneous(plane::Plane)
		# Convert a plane to homogeneous coordinates
		normal = plane.normal
		d = -dot(normal, plane.origin)
		return [normal; d]
	end

	function cylinder_rotation_from_axis(axis::Vector{<:Number})
		# Create a rotation matrix that aligns the Z-axis with the given axis
		# Assuming axis is a unit vector
		axis = normalize(axis)
		z_axis = [0.0, 0.0, 1.0]
		rotation = rotation_between_lines(z_axis, axis)
		return rotation
	end

	function get_tangentpoints_circle_point(circle::Circle, point::Vector{<:Number})
		variation = point - circle.center
		orthogonal_variation = variation
		axis = normalize(circle.axis)
		if !isnothing(axis)
			orthogonal_variation = cross(axis, variation)
		else
			orthogonal_variation = variation .* [1, -1]
		end
		d = norm(variation)
		if d >= circle.radius
			rho = circle.radius/d
			ad = rho^2
			bd = rho*sqrt(1-rho^2)
			T1 = circle.center + ad * variation + bd * orthogonal_variation
			T2 = circle.center + ad * variation - bd * orthogonal_variation

			if (d/circle.radius-1) ≃ 1E-8
				throw(TangentLineNotFound("The point is on the circle"))
			else
				return (T1, T2)
			end
		else
			throw(TangentLineNotFound("The point is inside the circle"))
		end

		throw(TangentLineNotFound("No tangent line possible"))
	end

	function get_cylinder_contours(cylinder::CylinderProperties, camera::CameraProperties)
		P = camera.matrix
		iCylinder = cylinder.dual_matrix
		iPlane = inv(cylinder.transform) * Diagonal([0, 0, 1, 0]) * inv(cylinder.transform')
		# Project the 3D dual quadrics to 2D
		projected_cylinder = P * iCylinder * transpose(P)
		projectedPlane = P * iPlane * transpose(P)
		
		# Get basis vectors for the line at infinity
		_, _, V = svd(projectedPlane)
		baseVec1 = V[:, 2]
		baseVec2 = V[:, 3]

		# Build the quadratic equation coefficients
		coeffs = [
			dot(baseVec1, projected_cylinder * baseVec1),
			2 * dot(baseVec1, projected_cylinder * baseVec2),
			dot(baseVec2, projected_cylinder * baseVec2)
		]

		# Solve the quadratic equation
		p = Polynomial(reverse(coeffs))  # Reverse because Polynomial expects [c, b, a]
		rootsOfConic = roots(p)
		
		# Compute the two tangent lines
		tangentLine1 = baseVec1 * rootsOfConic[1] + baseVec2
		tangentLine2 = baseVec1 * rootsOfConic[2] + baseVec2
	
		return tangentLine1, tangentLine2
	end
end