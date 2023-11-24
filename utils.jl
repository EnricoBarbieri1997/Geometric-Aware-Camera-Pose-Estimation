module Utils
	export randRange

	using Random

	function randRange(a, b)
		return rand(Float64) * (b - a) + a
	end
	
	function randRange(range::Tuple{Number, Number})
		return randRange(range...)
	end
end