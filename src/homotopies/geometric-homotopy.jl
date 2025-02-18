export GeometricHomotopy

using ..Geometry: homogeneous_anglebetween
using LinearAlgebra: cross
using HomotopyContinuation

"""
    GeometricHomotopy(F::Union{AbstractSystem,System}; start_parameters, target_parameters)
    GeometricHomotopy(F::Union{AbstractSystem,System}, start_parameters, target_parameters)

Construct the parameter homotopy ``H(x,t) = F(x; t p + (1 - t) q)`` where ``p`` is
`start_parameters` and ``q`` is `target_parameters`.
"""
struct GeometricHomotopy{T<:AbstractSystem} <: AbstractHomotopy
    F::T
    p::Vector{ComplexF64}
    q::Vector{ComplexF64}
    #cache
	angles_between::Vector{Real}
    intersections::Vector{Vector{Real}}
    t_cache::Base.RefValue{ComplexF64}
    pt::Vector{ComplexF64}
    taylor_pt::TaylorVector{2,ComplexF64}
end

function GeometricHomotopy(
    F;
    start_parameters::AbstractVector,
    target_parameters::AbstractVector,
)
    GeometricHomotopy(F, start_parameters, target_parameters)
end
function GeometricHomotopy(
    F::ModelKit.System,
    p::AbstractVector,
    q::AbstractVector;
    compile::Union{Bool,Symbol} = true,
)
    GeometricHomotopy(fixed(F; compile = compile), p, q)
end
function GeometricHomotopy(F::AbstractSystem, p::AbstractVector, q::AbstractVector)
    @assert length(p) == length(q) == nparameters(F)

    p̂ = Vector{ComplexF64}(p)
    q̂ = Vector{ComplexF64}(q)
    taylor_pt = TaylorVector{2}(ComplexF64, length(q))
    pt = copy(p̂)

    number_of_lines = convert(Int, length(p) / 3)

    angles_between = zeros(Real, number_of_lines)
    intersections = Vector{Vector{Real}}(undef, number_of_lines)
    for i in 1:number_of_lines
        index = (i-1)*3+1
        line_start = p[index:index+2]
        line_start = line_start / line_start[3]
        line_target = q[index:index+2]
        line_target = line_target / line_target[3]
        angles_between[i] = homogeneous_anglebetween(line_start, line_target)
        intersections[i] = cross(line_start, line_target)
    end

    GeometricHomotopy(F, p̂, q̂, angles_between, intersections, Ref(complex(NaN)), pt, taylor_pt)
end

Base.size(H::GeometricHomotopy) = size(H.F)

function start_parameters!(H::GeometricHomotopy, p)
    H.p .= p
    # void cache
    H.t_cache[] = NaN
    H
end
function target_parameters!(H::GeometricHomotopy, q)
    H.q .= q
    H.t_cache[] = NaN
    H
end
function parameters!(H::GeometricHomotopy, p, q)
    H.p .= p
    H.q .= q
    H.t_cache[] = NaN
    H
end

function tp!(H::GeometricHomotopy, tinput::Union{ComplexF64,Float64})
    tinput == H.t_cache[] && return H.taylor_pt
    t = if imag(tinput) == 0 real(tinput) else tinput end

    parameters = zeros(convert(Int, length(H.p)))
    for i in 1:length(H.angles_between)
        index = (i-1)*3+1
        line_start = H.p[index:index+2]
        angle = H.angles_between[i] * (1.0-t)
        intersection = (H.intersections[i] / H.intersections[i][3])[1:2]
        rotation = [
            cos(angle) -sin(angle);
            sin(angle) cos(angle);
        ]
        transform = vcat(hcat(rotation, -rotation * intersection + intersection), [0, 0, 1]')
        parameters[index:index+2] = transform * line_start
    end

    @inbounds for i = 1:length(H.taylor_pt)
        ptᵢ = parameters[i]
        H.pt[i] = ptᵢ
        H.taylor_pt[i] = (ptᵢ, H.p[i] - H.q[i])
    end
    H.t_cache[] = t

    H.taylor_pt
end

function ModelKit.evaluate!(u, H::GeometricHomotopy, x, t)
    tp!(H, t)
    a = evaluate!(u, H.F, x, H.pt)
    display(a)
    a
end

function ModelKit.evaluate_and_jacobian!(u, U, H::GeometricHomotopy, x, t)
    tp!(H, t)
    evaluate_and_jacobian!(u, U, H.F, x, H.pt)
end

function ModelKit.taylor!(u, v::Val, H::GeometricHomotopy, tx, t)
    taylor!(u, v, H.F, tx, tp!(H, t))
    u
end