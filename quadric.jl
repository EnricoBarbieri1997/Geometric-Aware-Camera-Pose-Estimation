module Quadric
	export ToFormula

	indexToVariable = Dict{Int, String}()
	indexToVariable[1] = "x"
	indexToVariable[2] = "y"
	indexToVariable[3] = "z"
	indexToVariable[4] = ""

	function ToFormula(quadric::Matrix)
		formula = ""
		for i in 1:4
			for j in i:4
				coefficient = i == j ? quadric[i, j] : (2 * quadric[i, j])
				variable1 = indexToVariable[i]
				variable2 = indexToVariable[j]
				variables = i == j ? (variable1 == "" ? variable1 : (variable1 * "^2")) : (variable1 * variable2)
				formula = formula * "(" * string(coefficient) * ")" * (variables == "" ? variables : ("*" * variables )) * " + "
			end
		end
		formula = formula[1:end-3] * " = 0"
		return formula
	end
end