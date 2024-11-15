module Geometry
	export Plane, Line, Point, Circle, Cylinder, issame_line, project_point_into_line, project_point_into_plane, get_tangentpoints_circle_point, get_cylinder_contours

	using ..Utils
	using LinearAlgebra: cross, dot, norm, normalize

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
	Cylinder = Circle

	function issame_line(line₁::Line, line₂::Line)
		direction_factors = line₁.direction ./ line₂.direction
		issame_direction = allequal(direction_factors)
		if !issame_direction
			return false
		end

		origin_factors = (line₁.origin - line₂.origin) ./ line₁.direction
		return allequal(origin_factors)
	end

	function project_point_into_line(point::Vector{<:Number}, line::Line)
		direction = line.direction
		origin = line.origin
		v = point - origin
		return origin + dot(v, direction) / dot(direction, direction) * direction
	end

	function project_point_into_plane(point::Vector{<:Number}, plane::Plane)
		normal = plane.normal
		origin = plane.origin
		v = point - origin
		d = dot(v, normal)
		return point - d * normal
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
				throw(ArgumentError("The point is on the circle"))
			else
				return (T1, T2)
			end
		else
			throw(ArgumentError("The point is inside the circle"))
		end

		throw(ArgumentError("No tangent line possible"))
	end

	function get_cylinder_contours(cylinder::Cylinder, point::Vector{<:Number}, cameraMatrix::Matrix{<:Number})
		circlecenter = project_point_into_line(point, Line(cylinder.center, cylinder.axis))
		tangentpoint₁, tangentpoint₂ = get_tangentpoints_circle_point(
			Circle(circlecenter, cylinder.radius, cylinder.axis),
			point
		)
		projected_tangentpoint₁ = cameraMatrix * [tangentpoint₁; 1]
		projected_tangentpoint₂ = cameraMatrix * [tangentpoint₂; 1]
		projected_cylinderaxis = cameraMatrix * [cylinder.axis; 0]
		contour₁ = Line(projected_tangentpoint₁, projected_cylinderaxis - projected_tangentpoint₁)
		contour₂ = Line(projected_tangentpoint₂, projected_cylinderaxis - projected_tangentpoint₂)

		return (contour₁, contour₂)
	end
end