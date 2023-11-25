module Utils
	export randRange

	using Random

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