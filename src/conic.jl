module Conic
	export ConicProperties

	struct ConicProperties
		matrix::Matrix{<:Number}
		singular_point::Vector{Number}
		dual_matrix::Matrix{<:Number}
	end
end