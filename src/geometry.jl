module Geometry
	export Plane, Line, Point, Circle, Cylinder, TangentLineNotFound, issame_line, homogeneous_to_line, line_to_homogenous, homogeneous_line_intercept, homogeneous_anglebetween, project_point_into_line, project_point_into_plane, get_tangentpoints_circle_point, get_cylinder_contours

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

	struct TangentLineNotFound <: Exception
		msg::String
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
				throw(TangentLineNotFound("The point is on the circle"))
			else
				return (T1, T2)
			end
		else
			throw(TangentLineNotFound("The point is inside the circle"))
		end

		throw(TangentLineNotFound("No tangent line possible"))
	end

	function get_cylinder_contours(cylinder::Cylinder, cameraCenter::Vector{<:Number}, cameraMatrix::Matrix{<:Number})
		circlecenter = project_point_into_line(cameraCenter, Line(cylinder.center, cylinder.axis))
		tangentpoint₁, tangentpoint₂ = get_tangentpoints_circle_point(
			Circle(circlecenter, cylinder.radius, cylinder.axis),
			cameraCenter
		)
		projected_tangentpoint₁ = cameraMatrix * [tangentpoint₁; 1]
		projected_tangentpoint₂ = cameraMatrix * [tangentpoint₂; 1]
		projected_cylinderaxis = cameraMatrix * [cylinder.axis; 0]

		projected_tangentpoint₁ = projected_tangentpoint₁ ./ projected_tangentpoint₁[3]
		projected_tangentpoint₂ = projected_tangentpoint₂ ./ projected_tangentpoint₂[3]
		projected_cylinderaxis = projected_cylinderaxis ./ projected_cylinderaxis[3]

		projected_tangentpoint₁ = projected_tangentpoint₁[1:2]
		projected_tangentpoint₂ = projected_tangentpoint₂[1:2]
		projected_cylinderaxis = projected_cylinderaxis[1:2]

		contour₁ = Line(projected_cylinderaxis, projected_tangentpoint₁ - projected_cylinderaxis)
		contour₂ = Line(projected_cylinderaxis, projected_tangentpoint₂ - projected_cylinderaxis)

		return (contour₁, contour₂)
	end
end