module Utils
	export almostEqual, ≃, randRange

	using Random

	function almostEqual(x::Number, y::Number)
		return abs(x - y) < 1e-6
	end

	function almostEqual(x::Matrix, y::Number)
		if (size(x) == (1,1))
			return almostEqual(x[1][1], y)
		end
		return false
	end

	function almostEqual(x::Number, y::Matrix)
		if (size(y) == (1,1))
			return almostEqual(x, y[1][1])
		end
		return false
	end

	function almostEqual(x::Matrix, y::Matrix)
		if (size(x) == (1,1) && size(y) == (1,1))
			return almostEqual(x[1][1], y[1][1])
		end
		if (!(size(x) == size(y)))
			return false
		end
		for i in eachindex(x)
			if (!almostEqual(x[i], y[i]))
				return false
			end
		end
		return true
	end

	function almostEqual(x::Vector, y::Number)
		if (size(x) == (1,))
			return almostEqual(x[1], y)
		end
		return false
	end

	function almostEqual(x::Number, y::Vector)
		if (size(y) == (1,))
			return almostEqual(x, y[1])
		end
		return false
	end

	function almostEqual(x::Vector, y::Vector)
		if length(x) != length(y)
			return false
		end
		for i in eachindex(x)
			if !almostEqual(x[i], y[i])
				return false
			end
		end
		return true
	end
	
	function ≃(x::Number, y::Number)
		return almostEqual(x, y)
	end

	function ≃(x::Matrix, y::Number)
		return almostEqual(x, y)
	end

	function ≃(x::Number, y::Matrix)
		return almostEqual(x, y)
	end

	function ≃(x::Matrix, y::Matrix)
		return almostEqual(x, y)
	end

	function ≃(x::Vector, y::Number)
		return almostEqual(x, y)
	end

	function ≃(x::Number, y::Vector)
		return almostEqual(x, y)
	end

	function ≃(x::Vector, y::Vector)
		return almostEqual(x, y)
	end

	function randRange(a::Number, b::Number):: Float64
		return rand(Float64) * (b - a) + a
	end
	
	function randRange(range::Tuple{Number, Number}):: Float64
		return randRange(range...)
	end

	function randRange(a, b, n):: Array{Float64}
		rands = []
		for i in 1:n
			push!(rands, randRange(a, b))
		end
		return rands
	end

	function randRange(range::Tuple{Number, Number}, n):: Array{Float64}
		return randRange(range..., n)
	end

	function randRange(range::Vector{Tuple{Float64, Float64}}):: Array{Float64}
		rands = []
		for r in range
			push!(rands, randRange(r))
		end
		return rands
	end

	function randRange(range::Vector{Tuple{Int64, Int64}}):: Array{Float64}
		return randRange(map(r -> (Float64(r[1]), Float64(r[2])), range))
	end
end