module Geometry
	export get_tangentline_circle_point

	using ..Utils

	struct Point
		x::Number
		y::Number
	end
	struct Circle
		center::Vecotr{Number}
		radius::Number
		axis::Union{Vecotr{Number}, Nothing}
	end
	struct Plane
		origin::Vector{Number}
		normal::Vector{Number}
	end

	function get_tangentpoints_circle_point(circle::Circle, point::Vector{Number})
		variation = point - circle.center
		orthogonal_variation = variation
		if !isnothing(circle.axis)
			orthogonal_variation = cross(circle.axis, variation)
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
				throw(ArgumentError("The point is on the circle"))
			else
				return (T1, T2)
			end
		else
			throw(ArgumentError("The point is inside the circle"))
		end

		throw(ArgumentError("No tangent line possible"))
	end

	function project_point_into_plane(point::Vector{Number}, plane::Plane)
		normal = plane.normal
		origin = plane.origin
		v = point - origin
		d = dot(v, normal)
		return point - d * normal
	end
end