module Formulas
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

	module Conic
		export ToFormula

		indexToVariable = Dict{Int, String}()
		indexToVariable[1] = "x"
		indexToVariable[2] = "y"
		indexToVariable[3] = ""

		function ToFormula(conic::Matrix)
			formula = ""
			for i in 1:3
				for j in i:3
					coefficient = i == j ? conic[i, j] : (2 * conic[i, j])
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

	module Plane 
		export ToFormula

		function ToFormula(plane::Array)
			return string(plane[1]) * "x + " * string(plane[2]) * "y + " * string(plane[3]) * "z + " * string(plane[4]) * " = 0"
		end
	end

	module PointFormulas
		export ToFormula

		function ToFormula(point::Vector)
			return "(" * string(point[1]) * ", " * string(point[2]) * ", " * string(point[3]) * ")"
		end
	end
end