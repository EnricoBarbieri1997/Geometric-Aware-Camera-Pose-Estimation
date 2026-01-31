module EquationSystems
export stack_homotopy_parameters, build_intrinsic_rotation_conic_system, build_intrinsic_rotation_translation_conic_system, build_intrinsic_rotation_translation_conic_system_calibrated, variables_jacobian_rank, joint_jacobian_rank

using ..CylindersBasedCameraResectioning: IMAGE_HEIGHT, IMAGE_WIDTH
using ..Camera: IntrinsicParameters, build_camera_matrix
using ..Space: build_rotation_matrix

using HomotopyContinuation
using LinearAlgebra: det, I, rank

module Problems
using ....Camera: CameraProperties
module IntrinsicParameters
@enum T begin
    focal_length_x = 0b00001
    focal_length_y = 0b00010
    skew = 0b00100
    principal_point_x = 0b01000
    principal_point_y = 0b10000
end
function Base.:|(a::T, b::T)
    UInt8(a) | UInt8(b)
end
function Base.:&(a::T, b::T)
    UInt8(a) & UInt8(b)
end
function Base.:|(a::Any, b::T)
    a | UInt8(b)
end
function Base.:&(a::Any, b::T)
    a & UInt8(b)
end
function Base.:|(a::T, b::Any)
    UInt8(a) | b
end
function Base.:&(a::T, b::Any)
    UInt8(a) & b
end
module has
using ..IntrinsicParameters: focal_length_x, focal_length_y, skew as skewParameter, principal_point_x, principal_point_y
function fₓ(config::Any)
    return (UInt8(config) & focal_length_x) != 0
end
function fᵧ(config::Any)
    return (UInt8(config) & focal_length_y) != 0
end
function skew(config::Any)
    return (UInt8(config) & skewParameter) != 0
end
function cₓ(config::Any)
    return (UInt8(config) & principal_point_x) != 0
end
function cᵧ(config::Any)
    return (UInt8(config) & principal_point_y) != 0
end
end
module Configurations
using ..IntrinsicParameters: focal_length_x, focal_length_y, skew as skew_val, principal_point_x, principal_point_y
@enum T begin
    none = 0
    fₓ = UInt8(focal_length_x)
    fᵧ = UInt8(focal_length_y)
    skew = UInt8(skew_val)
    fₓ_fᵧ = focal_length_x | focal_length_y
    cₓ_cᵧ = principal_point_x | principal_point_y
    fₓ_fᵧ_skew = focal_length_x | focal_length_y | skew_val
    fₓ_fᵧ_skew_cₓ = focal_length_x | focal_length_y | skew_val | principal_point_x
    fₓ_fᵧ_skew_cᵧ = focal_length_x | focal_length_y | skew_val | principal_point_y
    fₓ_fᵧ_cₓ_cᵧ = focal_length_x | focal_length_y | principal_point_x | principal_point_y
    fₓ_fᵧ_skew_cₓ_cᵧ = focal_length_x | focal_length_y | skew_val | principal_point_x | principal_point_y
end
end
end
mutable struct CylinderCameraContoursProblemValidationData
    lines::Array{Float64,2}
    points_at_infinity::Array{Float64,2}
    dualquadrics::Array{Float64,3}
    line_indexes::Array{Float64,1}
end
mutable struct CylinderCameraContoursProblem
    camera::CameraProperties
    lines::Array{Float64,2}
    noise_free_lines::Array{Float64,2}
    points_at_infinity::Array{Float64,2}
    dualquadrics::Array{Float64,3}
    line_indexes::Array{Float64,1}
    validation::CylinderCameraContoursProblemValidationData
    intrinsic_configuration::UInt8
end

function CylinderCameraContoursProblem(
    camera::CameraProperties,
    lines::Array{Float64,2},
    noise_free_lines::Array{Float64,2},
    points_at_infinity::Array{Float64,2},
    dualquadrics::Array{Float64,3},
    intrinsic_configuration::UInt8=IntrinsicParameters.focal_length_x | IntrinsicParameters.focal_length_y | IntrinsicParameters.skew | IntrinsicParameters.principal_point_x | IntrinsicParameters.principal_point_y
)
    println(typeof(camera))
    println(typeof(lines))
    println(typeof(noise_free_lines))
    println(typeof(points_at_infinity))
    println(typeof(dualquadrics))
    println(typeof(intrinsic_configuration))

    return CylinderCameraContoursProblem(
        camera,
        lines,
        noise_free_lines,
        points_at_infinity,
        dualquadrics,
        collect(1:size(lines)[1]),
        CylinderCameraContoursProblemValidationData(
            [], [], [], []
        ),
        intrinsic_configuration
    )
end
end

function stack_homotopy_parameters(parameters...)
    stacked_parameters = []
    for parameter in parameters
        stacked_parameters = vcat(stacked_parameters, vec(parameter))
    end
    return stacked_parameters
end

function add_rotation_constraints!(system_to_solve, R)
    push!(system_to_solve, det(R) - 1)
    for eq in vec(R * R' - I)
        push!(system_to_solve, eq)
    end
end

function build_intrinsic_rotation_conic_system(
    problems::Vector{Problems.CylinderCameraContoursProblem};
)
    problems_count = length(problems)
    if problems_count == 0
        throw(ArgumentError("At least one problem is needed"))
    end

    default_intrinsic = problems[1].camera.intrinsic
    factor = 1.0 / 3000.0
    fₓ = default_intrinsic[1, 1] * factor
    fᵧ = default_intrinsic[2, 2] * factor
    skew = default_intrinsic[1, 2] * factor
    cₓ = default_intrinsic[1, 3] * factor
    cᵧ = default_intrinsic[2, 3] * factor

    system_to_solve = []
    variables::Vector{HomotopyContinuation.ModelKit.Variable} = []
    parameters::Vector{HomotopyContinuation.ModelKit.Variable} = []

    intrinsic_configuration = problems[1].intrinsic_configuration
    if Problems.IntrinsicParameters.has.fₓ(intrinsic_configuration)
        @var fₓ
        push!(variables, fₓ)
    end
    if Problems.IntrinsicParameters.has.fᵧ(intrinsic_configuration)
        @var fᵧ
        push!(variables, fᵧ)
    end
    if Problems.IntrinsicParameters.has.skew(intrinsic_configuration)
        @var skew
        push!(variables, skew)
    end
    if Problems.IntrinsicParameters.has.cₓ(intrinsic_configuration)
        @var cₓ
        push!(variables, cₓ)
    end
    if Problems.IntrinsicParameters.has.cᵧ(intrinsic_configuration)
        @var cᵧ
        push!(variables, cᵧ)
    end

    intrinsic = [
        fₓ skew cₓ;
        0 fᵧ cᵧ;
        0 0 factor
    ]

    intrinsic_count = count_ones(UInt8(intrinsic_configuration))
    extrinsic_count = size(problems)[1] * 3
    lines_count = intrinsic_count + extrinsic_count
    lines_per_problem = zeros(Int64, problems_count)

    problem_index = 1
    for i in 1:lines_count
        lines_per_problem[problem_index] += 1
        problem_index += 1
        if problem_index > problems_count
            problem_index = 1
        end
    end

    for (index, problem) in enumerate(problems)
        lines_to_pick = lines_per_problem[index]
        Rparams = [
            Variable("R$(index)", i) for i in 1:3
        ]
        R = build_rotation_matrix(Rparams..., false)
        lines = reshape([
                Variable("lines$(index)", i, j)
                for i in 1:lines_to_pick, j in 1:3
            ], lines_to_pick, 3)
        display(lines)
        for line_index in 1:lines_to_pick
            # if index == 2 && line_index == 2
            #     display("params")
            #     display(lines[line_index, :])
            #     display(intrinsic)
            #     display(R)
            #     display(problem.points_at_infinity[line_index, :])
            # end
            equation = lines[line_index, :]' * intrinsic * R * problem.points_at_infinity[line_index, :]
            push!(system_to_solve, equation)
        end

        variables = stack_homotopy_parameters(variables, Rparams)
        parameters = stack_homotopy_parameters(parameters, lines)
        display(parameters)
    end

    a = System(system_to_solve, variables=variables, parameters=parameters)
    return a
end

function build_intrinsic_rotation_translation_conic_system(problem::Problems.CylinderCameraContoursProblem)
    lines_count = 3
    @var tx ty tz
    @var lines[1:lines_count, 1:3]
    P = build_camera_matrix(
        problem.camera.intrinsic ./ problem.camera.intrinsic[2, 2],
        problem.camera.quaternion_rotation,
        [tx, ty, tz]
    )

    system_to_solve = []
    for i in 1:lines_count
        equation = (lines[i, :]' * P * problem.dualquadrics[i, :, :] * P' * lines[i, :]) / (IMAGE_HEIGHT * IMAGE_WIDTH)
        push!(system_to_solve, equation)
    end
    parameters::Vector{HomotopyContinuation.ModelKit.Variable} = stack_homotopy_parameters(lines)
    return System(system_to_solve, variables=[tx, ty, tz], parameters=parameters)
end

function build_intrinsic_rotation_translation_conic_system_calibrated(problem::Problems.CylinderCameraContoursProblem)
    new_problem = deepcopy(problem)
    new_problem.camera.intrinsic = Matrix{Float64}(I, 3, 3)
    return build_intrinsic_rotation_translation_conic_system(new_problem)
end

function variables_jacobian(F::System, solution, parameters)
    jacobian = differentiate(F.expressions, F.variables)
    return evaluate(jacobian, F.variables => solution, F.parameters => parameters)
end

function variables_jacobian_rank(F::System, solution, parameters)
    return rank(variables_jacobian(F, solution, parameters))
end

function parameters_jacobian(F::System, solution, parameters)
    jacobian = differentiate(F.expressions, F.parameters)
    return evaluate(jacobian, F.variables => solution, F.parameters => parameters)
end

function parameters_jacobian_rank(F::System, solution, parameters)
    return rank(parameters_jacobian(F, solution, parameters))
end

function joint_jacobian_rank(F::System, solution, parameters)
    return rank(hcat(parameters_jacobian(F, solution, parameters), variables_jacobian(F, solution, parameters)))
end

module Minimization
using ..EquationSystems
using ..Problems

using HomotopyContinuation
using LinearAlgebra: norm

function bestsolution(system, solutions, parameters)
    best_solution = nothing
    best_residual = Inf
    for solution in solutions
        residual = norm(system(solution, parameters))
        if residual < best_residual
            best_residual = residual
            best_solution = solution
        end
    end
    return best_solution
end
function build_intrinsic_rotation_conic_system(
    problems::Vector{Problems.CylinderCameraContoursProblem};
    args...,
)
    sys = EquationSystems.build_intrinsic_rotation_conic_system(problems; args...)
    expressions = sys.expressions
    minimizer = sum(expressions .^ 2)
    variables = sys.variables
    parameters = sys.parameters

    diff = expand.(differentiate(minimizer, variables))

    return System(diff, variables=variables, parameters=parameters)
end
end

module SingleProblem
function build_intrinsic_rotation_conic_system(lines_values::Matrix{<:Number})
end
function build_intrinsic_rotation_translation_conic_system(intrinsic, rotation, lines_values)
end
end
end
