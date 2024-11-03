module CylindersBasedCameraResectioning
    include("includes.jl")

    using .Space: Transformation, RandomTransformation, IdentityTransformation, build_rotation_matrix
    using .Camera: CameraMatrix
    using .Plotting: initFigure, plot2DPoints, Plot3DCameraInput, plot3DCamera, Plot3DCylindersInput, plot3DCylinders, plot2DCylinders
    using .Debug
    using .Utils
    using LinearAlgebra: deg2rad, diagm, dot, normalize, svd, pinv
    using HomotopyContinuation, Polynomials, Rotations
    using Random

    function main()
        Random.seed!(7)

        numberOfCylinders = 4
        cylinders = Array{Tuple{Matrix{Float64}, Tuple{Matrix{Float64}, Vector{Float64}}}}(undef, numberOfCylinders)
        transforms = Array{Matrix{Float64}}(undef, numberOfCylinders)
        radiuses = Array{Tuple{Number, Number}}(undef, numberOfCylinders)
        dualSingularPlanes = Array{Matrix{Float64}}(undef, numberOfCylinders)
        for i in 1:numberOfCylinders
            transforms[i] = RandomTransformation()
            # transforms[i] = IdentityTransformation()
            radius = randRange((1, 3), 2)
            # radius = (1, 1)
            radiuses[i] = (radius[1], radius[2])
            cylinders[i] = Cylinder.StandardAndDual(transforms[i], radiuses[i])

            begin #asserts
                # @assert cylinders[i][1] * cylinders[i][2][1] ≃ diagm([1, 1, 0, 1]) "(-1) The dual quadric is correct"

                @assert cylinders[i][2][2]' * cylinders[i][1] * cylinders[i][2][2] ≃ 0 "(1) Singular point $(1) belongs to the cylinder $(1)"
                dualSingularPlanes[i] = inv(transforms[i]') * reshape([1, 0, 0, -radiuses[i][1]], :, 1)
                @assert (dualSingularPlanes[i]' * cylinders[i][2][1] * dualSingularPlanes[i]) ≃ 0 "(2) Perpendicular plane $(1) belongs to the dual cylinder $(1)"

                @assert (cylinders[i][1] * cylinders[i][2][2]) ≃ [0, 0, 0, 0] "(6) Singular point is right null space of cylinder matrix $(i)"

                @assert ((dualSingularPlanes[i]' * cylinders[i][2][2]) ≃ 0 && (dualSingularPlanes[i]' * cylinders[i][2][1] * dualSingularPlanes[i]) ≃ 0) "(7) Singular plane / point and dual quadric constraints $(i)"
                @assert cylinders[i][2][2][4] ≃ 0 "(10) Singular point is at infinity $(i)"
            end
        end

        cameraTranslation = (2.0, 30.0, 5.0)
        cameraRotation = (-83.0, 180.0, 0.0)
        # cameraTranslation = (0.0, 0.0, 0.0)
        # cameraRotation = (0.0, 0.0, 0.0)
        cameraPositionMatrix = Transformation(cameraTranslation, cameraRotation)
        cameraProjectionMatrix = CameraMatrix(cameraTranslation, cameraRotation, 2, 1)

        # display(cameraProjectionMatrix)
        # display(pinv(cameraProjectionMatrix))
        # display(pinv(cameraProjectionMatrix) * cameraProjectionMatrix)
        # display(cameraProjectionMatrix' * cameraProjectionMatrix)
        # display(cameraProjectionMatrix * cameraProjectionMatrix')

        conics = Array{Tuple{Matrix{Float64}, Tuple{Matrix{Float64}, Vector{Float64}}}}(undef, numberOfCylinders)
        for i in 1:numberOfCylinders
            conics[i] = (
                pinv(cameraProjectionMatrix') * cylinders[i][1] * pinv(cameraProjectionMatrix),
                (
                    cameraProjectionMatrix * cylinders[i][2][1] * cameraProjectionMatrix',
                    cameraProjectionMatrix * cylinders[i][2][2]
                )
            )

            begin #asserts
                # @assert conics[i][1] * conics[i][2][1] ≃ diagm([1, 1, 0]) "(-3) The dual conic is correct"
                cylinderProjection = cameraProjectionMatrix * cylinders[i][2][1] * cameraProjectionMatrix'
                @assert conics[i][2][1] ≃ cylinderProjection "(4) Dual conic $(1) is the transformation of the dual cylinder"
            end
        end

        function lines_from_conic(conic)
            @var x y z
            line = [x, y, z]
            conicSingularPoint = conic[2][2]
            conicMatrix = conic[2][1]
            f₁ = line' * conicSingularPoint
            f₂ = line' * conicMatrix * line
            f₃ = z - 1
            F = System([f₁, f₂, f₃])
            result = solve(F)
            lines = result
            return real_solutions(lines)
        end

        conicBorders = Array{Array{Vector{Float64}}}(undef, numberOfCylinders)
        for i in 1:numberOfCylinders
            lines = lines_from_conic(conics[i])
            conicBorders[i] = Array{Vector{Float64}}(undef, length(lines))
            for (j, line) in enumerate(lines)
                conicBorders[i][j] = line

                begin #asserts
                    @assert line' * conics[i][2][1] * line ≃ 0 "(3) Line of projected singular plane $(1) belongs to the dual conic $(1)"
                    @assert line' * cameraProjectionMatrix * cylinders[i][2][2] ≃ 0 "(8) Line $(j) of conic $(i) passes through the projected singular point"
                    @assert line' * cameraProjectionMatrix * cylinders[i][2][1] * cameraProjectionMatrix' * line ≃ 0 "(9) Line $(j) of conic $(i) is tangent to the projected cylinder"
                end
            end
        end

        singularPoints = Array{Tuple{Number, Number}}(undef, numberOfCylinders)
        for i in 1:numberOfCylinders
            singularPoint = conics[i][2][2]
            singularPoint = singularPoint ./ singularPoint[3]
            singularPoint = (singularPoint[1], singularPoint[2])
            singularPoints[i] = singularPoint
        end

        figure = initFigure()
        plot3DCamera(Plot3DCameraInput(
            cameraRotation,
            cameraTranslation
        ))
        plot3DCylinders(Plot3DCylindersInput(
            transforms,
            radiuses,
            numberOfCylinders,
            # cameraProjectionMatrix
        ))
        plot2DPoints(singularPoints)
        plot2DCylinders(conicBorders)

        # 3 line minimum to solve the pose
        numberOfLinesToSolveFor = 4
        lines = []
        pointAtInfinityToUse = []
        dualQuadricToUse = []
        possiblePicks = 1:(numberOfCylinders*2)
        pickedLines = []
        for _ in (1:numberOfLinesToSolveFor)
            lineIndex = rand(possiblePicks)
            possiblePicks = filter(x -> x != lineIndex, possiblePicks)
            i = ceil(Int, lineIndex / 2)
            j = (lineIndex - 1) % 2 + 1
            push!(pickedLines, (i, j))
            borders = conicBorders[i]
            line = borders[j]
            push!(lines, line)
            push!(pointAtInfinityToUse, cylinders[i][2][2][1:3])
            push!(dualQuadricToUse, cylinders[i][2][1])
        end

        @var x y z f
        Rₚ = build_rotation_matrix(x, y, z)

        systemToSolve = []
        for i in 1:numberOfLinesToSolveFor
            equation = lines[i]' * [
                f 0 0;
                0 f 0;
                0 0 1
            ] * Rₚ * pointAtInfinityToUse[i]
            push!(systemToSolve, equation)
        end

        F = System(systemToSolve, variables=[x, y, z, f])
        result = solve(F)
        display(result)

        # for solution in real_solutions(result)
        #     r = Rotations.QuatRotation(1, solution[1], solution[2], solution[3])
        #     for i in 1:numberOfLinesToSolveFor
        #         display((lines[i]' * r * pointAtInfinityToUse[i]) ≃ 0)
        #     end
        # end

        rotationCalculated = nothing
        focalLengthCalculated = nothing
        rotationSolutionError = Inf
        for (solutionIndex, solution) in enumerate(real_solutions(result))
            display("Trying solution $(solutionIndex)")
            xₛ = solution[1]
            yₛ = solution[2]
            zₛ = solution[3]
            fₛ = solution[4]
            rotation = Rotations.QuatRotation(1, xₛ , yₛ, zₛ)
            # rotation = buildRotationMatrix(xₛ, yₛ, zₛ, true)
            # display("Rotation angle: $(rotation_angle(rotation)), Rotation axis: $(rotation_axis(rotation))")
            acceptable = true
            currentError = 0
            for (i, border) in enumerate(conicBorders)
                for (j, line) in enumerate(border)
                    eq = line' * [
                        fₛ 0 0;
                        0 fₛ 0;
                        0 0 1
                    ] *rotation * cylinders[i][2][2][1:3]
                    currentError += abs(eq)
                    # display("Solution $(solutionIndex): $(eq)")
                    if (!(eq ≃ 0))
                        # display("Line $(j) of conic $(i) does not satisfy the constraints")
                        acceptable = false
                    end
                end
                if (!acceptable)
                    # break
                end
            end
            if (acceptable)
                # rotationCalculated = solution
                display("Solution $(solutionIndex) is acceptable")
                break
            end
            if (currentError < rotationSolutionError)
                rotationSolutionError = currentError
                rotationCalculated = rotation
                focalLengthCalculated = fₛ
            end
        end

        begin #asserts
            @assert rotationCalculated isa QuatRotation{Float64} "(11) Found rotation satisfies constraints"
        end

        # rot = deg2rad.(cameraRotation)
        # rotationCalculated = QuatRotation(RotXYZ(rot[1], rot[2], rot[3]))

        display("Error of the best rotation solution: $(rotationSolutionError)")
        display("Best solution for rotation: $(rad2deg.(Rotations.params(RotXYZ(rotationCalculated))))")
        display("Actual rotation: $(Rotations.params(RotXYZ(cameraRotation[1], cameraRotation[2], cameraRotation[3])))")

        function buildCameraMatrixForSolving(rotation, translation)
            return [focalLengthCalculated 0 0 0;
                0 focalLengthCalculated 0 0;
                0 0 focalLengthCalculated 0] * vcat(hcat(rotation, translation), [0 0 0 1])
        end

        @var tx ty tz
        P = buildCameraMatrixForSolving(rotationCalculated, [tx; ty; tz])

        systemToSolve = []
        for i in 1:3
            equation = lines[i]' * P * dualQuadricToUse[i] * P' * lines[i]
            push!(systemToSolve, equation)
        end

        F = System(systemToSolve, variables=[tx, ty, tz])
        result = solve(F)
        display(result)
        solutions = real_solutions(result)
        translationSolutionError = Inf
        translationCalculated = nothing
        for (solutionIndex, solution) in enumerate(solutions)
            txₛ = solution[1]
            tyₛ = solution[2]
            tzₛ = solution[3]
            translation = [txₛ, tyₛ, tzₛ]
            acceptable = true
            currentError = 0
            Pₛ = buildCameraMatrixForSolving(rotationCalculated, translation)
            for (i, conicBorder) in enumerate(conicBorders)
                for (j, line) in enumerate(conicBorder)
                    eq = line' * Pₛ * cylinders[i][2][1] * Pₛ' * line
                    currentError += abs(eq)
                    # display("Solution $(i): $(eq))")
                    # if ((i, j) in pickedLines)
                    #     display("Solution $(solutionIndex): $(eq ≃ 0)")
                    # end
                    if (!(eq ≃ 0))
                        # display("Line $(j) of conic $(i) does not satisfy the constraints")
                        acceptable = false
                    end
                    if (currentError < translationSolutionError)
                        translationSolutionError = currentError
                        translationCalculated = translation
                    end
                end
                if (!acceptable)
                    # break
                end
            end
        end

        display("Translation: $(translationCalculated)")

        calculatedCameraMatrix = buildCameraMatrixForSolving(rotationCalculated, translationCalculated)

        display("Camera projection matrix: $(cameraProjectionMatrix ./ cameraProjectionMatrix[3, 4])")
        display("Calculated projection camera matrix: $(calculatedCameraMatrix ./ calculatedCameraMatrix[3, 4])")

        reconstructedBorders = Array{Array{Vector{Float64}}}(undef, numberOfCylinders)
        for i in 1:numberOfCylinders
            reconstructedConic = (
                pinv(calculatedCameraMatrix') * cylinders[i][1] * pinv(calculatedCameraMatrix),
                (
                    calculatedCameraMatrix * cylinders[i][2][1] * calculatedCameraMatrix',
                    calculatedCameraMatrix * cylinders[i][2][2]
                )
            )
            local lines = lines_from_conic(reconstructedConic)
            reconstructedBorders[i] = Array{Vector{Float64}}(undef, length(lines))
            for (j, line) in enumerate(lines)
                reconstructedBorders[i][j] = line
            end
        end

        plot2DCylinders(reconstructedBorders; linestyle=:dash)
        figure
    end

    function generate_monodromy_solutions()
        
    end

    export main
end
