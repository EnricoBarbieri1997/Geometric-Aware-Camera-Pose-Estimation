export almostEqual, ≃

function almostEqual(x::Number, y::Number)
	return abs(x - y) < 1e-6
end

function ≃(x::Number, y::Number)
	return almostEqual(x, y)
end