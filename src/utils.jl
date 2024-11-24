module Utils
	export almostequal, ≃, rand_in_range, quat_from_rotmatrix, rotations_difference, translations_difference

	using LinearAlgebra: norm
	using Rotations
	using Random

	function almostequal(x::Number, y::Number)
		return abs(x - y) < 1e-6
	end

	function almostequal(x::Matrix, y::Number)
		if (size(x) == (1,1))
			return almostequal(x[1][1], y)
		end
		return false
	end

	function almostequal(x::Number, y::Matrix)
		if (size(y) == (1,1))
			return almostequal(x, y[1][1])
		end
		return false
	end

	function almostequal(x::Matrix, y::Matrix)
		if (size(x) == (1,1) && size(y) == (1,1))
			return almostequal(x[1][1], y[1][1])
		end
		if (!(size(x) == size(y)))
			return false
		end
		for i in eachindex(x)
			if (!almostequal(x[i], y[i]))
				return false
			end
		end
		return true
	end

	function almostequal(x::Vector, y::Number)
		if (size(x) == (1,))
			return almostequal(x[1], y)
		end
		return false
	end

	function almostequal(x::Number, y::Vector)
		if (size(y) == (1,))
			return almostequal(x, y[1])
		end
		return false
	end

	function almostequal(x::Vector, y::Vector)
		if length(x) != length(y)
			return false
		end
		for i in eachindex(x)
			if !almostequal(x[i], y[i])
				return false
			end
		end
		return true
	end
	
	function ≃(x::Number, y::Number)
		return almostequal(x, y)
	end

	function ≃(x::Matrix, y::Number)
		return almostequal(x, y)
	end

	function ≃(x::Number, y::Matrix)
		return almostequal(x, y)
	end

	function ≃(x::Matrix, y::Matrix)
		return almostequal(x, y)
	end

	function ≃(x::Vector, y::Number)
		return almostequal(x, y)
	end

	function ≃(x::Number, y::Vector)
		return almostequal(x, y)
	end

	function ≃(x::Vector, y::Vector)
		return almostequal(x, y)
	end

	function rand_in_range(a::Number, b::Number):: Float64
		return rand(Float64) * (b - a) + a
	end
	
	function rand_in_range(range::Tuple{Number, Number}):: Float64
		return rand_in_range(range...)
	end

	function rand_in_range(a, b, n):: Array{Float64}
		rands = []
		for i in 1:n
			push!(rands, rand_in_range(a, b))
		end
		return rands
	end

	function rand_in_range(range::Tuple{Number, Number}, n):: Array{Float64}
		return rand_in_range(range..., n)
	end

	function rand_in_range(range::Vector{Tuple{Float64, Float64}}):: Array{Float64}
		rands = []
		for r in range
			push!(rands, rand_in_range(r))
		end
		return rands
	end

	function rand_in_range(range::Vector{Tuple{Int64, Int64}}):: Array{Float64}
		return rand_in_range(map(r -> (Float64(r[1]), Float64(r[2])), range))
	end

	function best_solution(solutions::Vector{Float64}, tester::Function):: [Int64, Float64]
		best = 1
		currentBestError = Inf
		for (i, s) in enumerate(solutions)
			error = tester(s)
			if error < currentBestError
				best = i
				currentBestError = error
			end
		end
		return best
	end

	function quat_from_rotmatrix(dcm::AbstractMatrix{T}) where {T<:Real}
		a2 = 1 + dcm[1,1] + dcm[2,2] + dcm[3,3]
		a = sqrt(a2)/2
		b,c,d = (dcm[3,2]-dcm[2,3])/4a, (dcm[1,3]-dcm[3,1])/4a, (dcm[2,1]-dcm[1,2])/4a
		return QuatRotation(a,b,c,d)
	end

	function rotations_difference(q1::QuatRotation, q2::QuatRotation)
		return norm(Rotations.params(q1 * q2'))
	end

	function translations_difference(t1::Vector{<:Number}, t2::Vector{<:Number})
		return norm(t1 - t2)
	end
end