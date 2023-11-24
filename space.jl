module Space
	export Vec3, Vec3Tuple, Point, PointTuple

	struct Vec3
		x::Number
		y::Number
		z::Number
	end
	const Vec3Tuple = Tuple{Number, Number, Number}

	const Point = Vec3
	const PointTuple = Vec3Tuple
end
