module Utils
	export almostequal, ≃, rand_in_range, random_camera_lookingat_center, quat_from_rotmatrix, vector_difference, matrix_difference, rotations_difference, eulerangles_from_rotationmatrix, translations_difference, isvalid_startsolution

	using ..Camera: lookat_rotation

	using LinearAlgebra: diagm, norm, normalize, svdvals, tr
	using Rotations
	using Random
	using HomotopyContinuation: jacobian

	@enum EulerOrder begin
		EulerOrderXYZ = 1
		EulerOrderXZY = 2
		EulerOrderYXZ = 3
		EulerOrderYZX = 4
		EulerOrderZXY = 5
		EulerOrderZYX = 6
	end

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

	function vector_difference(v1::Vector{<:Number}, v2::Vector{<:Number})
		return norm(v1 - v2)
	end

	function normalized_diff(calculated, truth)
		if (calculated == 0 && truth == 0) return 0.0 end
		denominator = if (truth != 0) truth else calculated end
		return (abs(calculated - truth)/denominator)
	end

	function intrinsic_difference(calculated, truth)
		fₓ, fᵧ, cₓ, cᵧ, skew = calculated[1, 1], calculated[2, 2], calculated[1, 3], calculated[2, 3], calculated[1, 2]
		fₓₜ, fᵧₜ, cₓₜ, cᵧₜ, skewₜ = truth[1, 1], truth[2, 2], truth[1, 3], truth[2, 3], truth[1, 2]

		deltaF = normalized_diff(fₓ, fₓₜ)/2 + normalized_diff(fᵧ, fᵧₜ)/2
		deltaUV = normalized_diff(cₓ, cₓₜ)/2 + normalized_diff(cᵧ, cᵧₜ)/2
		deltaSkew = 2 * abs(skew - skewₜ)

		return [deltaF, deltaUV, deltaSkew]
	end

	function matrix_difference(m1, m2)
		return sqrt(sum((m1 - m2) .^ 2))
	end

	function rotations_difference(q1::QuatRotation, q2::QuatRotation)
		R1 = Matrix(q1)  # convert quaternion to rotation matrix
		R2 = Matrix(q2)
		diff = clamp((tr(R1 * transpose(R2)) - 1) / 2, -1, 1)
		return acosd(diff)
	end

	function translations_difference(t1::Vector{<:Number}, t2::Vector{<:Number})
		return norm(t1 - t2)
	end

	function isvalid_startsolution(system, solution, parameters)
		return minimum(svdvals(jacobian(system, solution, parameters))) > 1e-6
	end
	function eulerangles_from_rotationmatrix(rotation_matrix; order::EulerOrder = EulerOrderXYZ)
		r = rotation_matrix
		if (order == EulerOrderXYZ)
			sy = r[1,3]
			if (sy < 1)
				if (sy > -1)
					if (
						r[2,1] == 0 &&
						r[1,2] == 0 &&
						r[2,3] == 0 &&
						r[3,2] == 0 &&
						r[2,2] == 1
					)
						return [0, atan(r[1,3], r[1,1]), 0]
					else
						return [atan(-r[2,3], r[3,3]), asin(r[1,3]), atan(-r[1,2], r[1,1])]
					end
				else
					return [atan(r[2,1], r[2,2]), -π/2, 0]
				end
			else
				return [atan(r[2,1], r[2,2]), π/2, 0]
			end
		end
		if (order == EulerOrderXZY)
			sz = r[1,2]
			if (sz < 1)
				if (sz > -1)
					return [atan(r[3,2], r[2,2]), atan(r[1,3], r[1,1]), asin(-sz)]
				else
					# sz == -1
					return [-atan(r[2,3], r[3,3]), 0, π/2]
				end
			else
				# sz == 1
				return [-atan(r[2,3], r[3,3]), 0, -π/2]
			end
		end
		if (order == EulerOrderYXZ)
			m12 = r[2,3]
			if (m12 < 1)
				if (m12 > -1)
					if (
						r[2,1] == 0 &&
						r[1,2] == 0 &&
						r[1,3] == 0 &&
						r[3,1] == 0 &&
						r[1,1] == 1
					)
						return [atan(-m12, r[2,2]), 0, 0]
					else
						return [asin(-m12), atan(r[1,3], r[3,3]), atan(r[2,1], r[2,2])]
					end
				else
					return [π/2, atan(r[1,2], r[1,1]), 0]
				end
			else
				return [-π/2, -atan(r[1,2], r[1,1]), 0]
			end
		end
		if (order == EulerOrderYZX)
			sz = r[2,1]
			if (sz < 1)
				if (sz > -1)
					return [atan(-r[2,3], r[2,2]), atan(-r[3,1], r[1,1]), asin(sz)]
				else
					return [atan(r[3,2], r[3,3]), 0, -π/2]
				end
			else
				return [atan(r[3,2], r[3,3]), 0, π/2]
			end
		end
		if (order == EulerOrderZXY)
			sx = r[3,2]
			if (sx < 1)
				if (sx > -1)
					return [asin(sx), atan(-r[3,1], r[3,3]), atan(-r[1,2], r[2,2])]
				else
					return [-π/2, atan(r[1,3], r[1,1]), 0]
				end
			else
				return [π/2, atan(r[1,3], r[1,1]), 0]
			end
		end
		if (order == EulerOrderZYX)
			sy = r[3,1]
			if (sy < 1)
				if (sy > -1)
					return [atan(r[3,2], r[3,3]), asin(-sy), atan(r[2,1], r[1,1])]
				else
					return [0, π/2, -atan(r[1,2], r[2,2])]
				end
			else
				return [0, -π/2, -atan(r[1,2], r[2,2])]
			end
		end
		throw(ArgumentError("Invalid arguments"))
	end

	function random_camera_lookingat_center()
		camera_translationdirection = normalize(rand_in_range(-1.0, 1.0, 3))
		camera_translation = camera_translationdirection * rand_in_range(15.0, 30.0)
		camera_object_rotation = lookat_rotation(camera_translationdirection, [0, 0, 0], [0, 0, 1])
		return camera_translation, camera_object_rotation
	end
end