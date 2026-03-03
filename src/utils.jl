	module Utils
	export almostequal, ≃, rand_in_range, quat_from_rotmatrix, normalized_diff, vector_difference, matrix_difference, intrinsic_difference, rotations_difference, eulerangles_from_rotationmatrix, translations_difference, isvalid_startsolution, lines_clp_to_stack, get_view_CC_VP_P


	using LinearAlgebra: Adjoint, diagm, norm, normalize, svdvals, tr
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
		calculated_normalized = calculated / calculated[3,3]
		truth_normalized = truth / truth[3,3]
		fₓ, fᵧ, cₓ, cᵧ, skew = calculated_normalized[1, 1], calculated_normalized[2, 2], calculated_normalized[1, 3], calculated_normalized[2, 3], calculated_normalized[1, 2]
		fₓₜ, fᵧₜ, cₓₜ, cᵧₜ, skewₜ = truth_normalized[1, 1], truth_normalized[2, 2], truth_normalized[1, 3], truth_normalized[2, 3], truth_normalized[1, 2]

		deltaF = normalized_diff(fₓ, fₓₜ)/2 + normalized_diff(fᵧ, fᵧₜ)/2
		deltaUV = normalized_diff(cₓ, cₓₜ)/2 + normalized_diff(cᵧ, cᵧₜ)/2
		deltaSkew = abs(skew - skewₜ) # 2 * abs(skew - skewₜ)

		return [deltaF, deltaUV, deltaSkew]
	end

	function matrix_difference(m1, m2)
		return sqrt(sum((m1 - m2) .^ 2))
	end

	function rotations_difference(q1::Union{QuatRotation, Matrix{Float64}, Adjoint{Float64, Matrix{Float64}}}, q2::Union{QuatRotation, Matrix{Float64}, Adjoint{Float64, Matrix{Float64}}})
		R1 = Matrix(q1)  # convert quaternion to rotation matrix
		R2 = Matrix(q2)
		diff = clamp((tr(R1 * transpose(R2)) - 1) / 2, -1, 1)
		return acosd(diff)
	end

	function translations_difference(t1::Vector{<:Number}, t2::Vector{<:Number})
		return norm(t1 - t2) / max(norm(t1), norm(t2))
	end

	function isvalid_startsolution(system, solution, parameters)
		return minimum(svdvals(jacobian(system, solution, parameters))) > 1e-6
	end
	function eulerangles_from_rotationmatrix(rotation_matrix; order::EulerOrder = EulerOrderZYX)
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

	function get_view_CC_VP_P(CC, VP, P)
		cc = CC
		if length(cc) == 16
			cc = cc[[1, 2, 3, 4, 6, 7, 8, 11, 12, 16]]
		end
		if ndims(cc) == 1
			cc = reshape(cc, :, 1)
		end
		if size(cc, 1) != 10
			throw(ArgumentError("CC must be 10xN or a 4x4 matrix"))
		end

		vp = VP
		if ndims(vp) == 1
			vp = reshape(vp, :, 1)
		end
		if size(vp, 1) == 3
			vp = vcat(vp, zeros(eltype(vp), 1, size(vp, 2)))
		end
		if size(vp, 1) != 4
			throw(ArgumentError("VP must be 4xN or 3xN"))
		end

		vp = P * vp

		C1 = cc[1, :]
		C2 = cc[2, :]
		C3 = cc[3, :]
		C4 = cc[4, :]
		C5 = cc[5, :]
		C6 = cc[6, :]
		C7 = cc[7, :]
		C8 = cc[8, :]
		C9 = cc[9, :]
		C10 = cc[10, :]

		p1 = P[1]
		p2 = P[2]
		p3 = P[3]
		p4 = P[4]
		p5 = P[5]
		p6 = P[6]
		p7 = P[7]
		p8 = P[8]
		p9 = P[9]
		p10 = P[10]
		p11 = P[11]
		p12 = P[12]

		nc = size(cc, 2)
		conic = zeros(eltype(cc), 6, nc)
		conic[1, :] = p1 * (C1 * p1 + C2 * p4 + C3 * p7 + C4 * p10) + p4 * (C2 * p1 + C5 * p4 + C6 * p7 + C7 * p10) + p7 * (C3 * p1 + C6 * p4 + C8 * p7 + C9 * p10) + p10 * (C4 * p1 + C7 * p4 + C9 * p7 + C10 * p10)
		conic[2, :] = p1 * (C1 * p2 + C2 * p5 + C3 * p8 + C4 * p11) + p4 * (C2 * p2 + C5 * p5 + C6 * p8 + C7 * p11) + p7 * (C3 * p2 + C6 * p5 + C8 * p8 + C9 * p11) + p10 * (C4 * p2 + C7 * p5 + C9 * p8 + C10 * p11)
		conic[3, :] = p1 * (C1 * p3 + C2 * p6 + C3 * p9 + C4 * p12) + p4 * (C2 * p3 + C5 * p6 + C6 * p9 + C7 * p12) + p7 * (C3 * p3 + C6 * p6 + C8 * p9 + C9 * p12) + p10 * (C4 * p3 + C7 * p6 + C9 * p9 + C10 * p12)
		conic[4, :] = p2 * (C1 * p2 + C2 * p5 + C3 * p8 + C4 * p11) + p5 * (C2 * p2 + C5 * p5 + C6 * p8 + C7 * p11) + p8 * (C3 * p2 + C6 * p5 + C8 * p8 + C9 * p11) + p11 * (C4 * p2 + C7 * p5 + C9 * p8 + C10 * p11)
		conic[5, :] = p2 * (C1 * p3 + C2 * p6 + C3 * p9 + C4 * p12) + p5 * (C2 * p3 + C5 * p6 + C6 * p9 + C7 * p12) + p8 * (C3 * p3 + C6 * p6 + C8 * p9 + C9 * p12) + p11 * (C4 * p3 + C7 * p6 + C9 * p9 + C10 * p12)
		conic[6, :] = p3 * (C1 * p3 + C2 * p6 + C3 * p9 + C4 * p12) + p6 * (C2 * p3 + C5 * p6 + C6 * p9 + C7 * p12) + p9 * (C3 * p3 + C6 * p6 + C8 * p9 + C9 * p12) + p12 * (C4 * p3 + C7 * p6 + C9 * p9 + C10 * p12)

		return linesfromconic(conic, vp)
	end

	function linesfromconic(cc, vp)
		conic = cc
		if length(conic) == 9
			conic = conic[[1, 2, 3, 5, 6, 9]]
		end
		if ndims(conic) == 1
			conic = reshape(conic, :, 1)
		end
		if size(conic, 1) != 6
			throw(ArgumentError("Conic must be 6xN or a 3x3 matrix"))
		end

		vp_local = vp
		if ndims(vp_local) == 1
			vp_local = reshape(vp_local, :, 1)
		end
		if size(vp_local, 1) == 2
			vp_local = vcat(vp_local, ones(eltype(vp_local), 1, size(vp_local, 2)))
		end
		if size(vp_local, 1) != 3
			throw(ArgumentError("vp must be 3xN or 2xN"))
		end

		c1 = conic[1, :]
		c2 = conic[2, :]
		c3 = conic[3, :]
		c4 = conic[4, :]
		c5 = conic[5, :]
		c6 = conic[6, :]

		vpx = vp_local[1, :]
		vpy = vp_local[2, :]
		vpz = vp_local[3, :]

		C2 = c4 .* vpx.^ 2 - 2 * c2 .* vpx .* vpy + c1 .* vpy.^ 2
		C1 = 2 * c3 .* vpy.^ 2 - 2 * c5 .* vpx .* vpy + 2 * c4 .* vpx .* vpz - 2 * c2 .* vpy .* vpz
		C0 = c6 .* vpy.^ 2 - 2 * c5 .* vpy .* vpz + c4 .* vpz.^ 2

		discriminant = sqrt.(C1 .^ 2 - 4 * C0 .* C2)
		lx1 = -(C1 - discriminant) ./ (2 * C2)
		lx2 = -(C1 + discriminant) ./ (2 * C2)

		ly1 = -(vpz + lx1 .* vpx) ./ vpy
		ly2 = -(vpz + lx2 .* vpx) ./ vpy

		row = x -> reshape(x, 1, :)
		ll1 = vcat(row(lx1), row(ly1), ones(eltype(lx1), 1, size(conic, 2)))
		ll2 = vcat(row(lx2), row(ly2), ones(eltype(lx2), 1, size(conic, 2)))
		return ll1, ll2
	end
	# Cylinder_count x 2 x 3 to line_count x 3
	function lines_clp_to_stack(lines_clp::Array{Float64,3})::Array{Float64,2}
		cylinder_count = size(lines_clp, 1)
		line_count = cylinder_count * 2
		return reshape(permutedims(lines_clp, (2, 1, 3)), line_count, 3)
	end
end
