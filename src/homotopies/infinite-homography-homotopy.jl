export InfiniteHomographyHomotopy, interpolate_line, verify_line_through_vanishing_point, matrix_exp, matrix_log, tp!

using LinearAlgebra: cross, norm, normalize, dot, pinv, I, eigen, Diagonal
using HomotopyContinuation

"""
    pencil_basis(v)

Compute two independent lines e, f that pass through point v (i.e., eᵀv = 0 and fᵀv = 0).
These form a basis for the pencil of lines through v.
"""
function pencil_basis(v::AbstractVector)
    v = normalize(v)
    vx, vy, vw = v[1], v[2], v[3]

    # First basis line: perpendicular in x-y plane
    if abs(vx) > abs(vy)
        e = normalize([-vy, vx, 0.0])
    else
        e = normalize([vy, -vx, 0.0])
    end

    # Second basis line: orthogonal to both v and e in line space
    # Use cross product in the dual space
    f = cross(v, e)
    f = normalize(f)

    # Verify: eᵀv ≈ 0 and fᵀv ≈ 0
    @assert abs(dot(e, v)) < 1e-10 "e must pass through v"
    @assert abs(dot(f, v)) < 1e-10 "f must pass through v"

    return e, f
end

"""
    express_in_pencil_basis(line, e, f)

Express a line in the pencil basis (e, f).
Returns (λ, μ) such that line ≈ λ*e + μ*f.
"""
function express_in_pencil_basis(line::AbstractVector, e::AbstractVector, f::AbstractVector)
    # Solve [e f] * [λ; μ] = line via pseudo-inverse
    basis_matrix = hcat(e, f)  # 3x2
    coeffs = pinv(basis_matrix) * line  # 2x1
    return coeffs[1], coeffs[2]
end

"""
    matrix_log(H)

Compute the matrix logarithm of a 3x3 matrix H.
Uses eigendecomposition: if H = V * D * V⁻¹, then log(H) = V * log(D) * V⁻¹.
"""
function matrix_log(H::AbstractMatrix)
    F = eigen(H)
    V = F.vectors
    D = Diagonal(log.(Complex.(F.values)))
    return real.(V * D * inv(V))
end

"""
    matrix_exp(A)

Compute the matrix exponential of a 3x3 matrix A.
"""
function matrix_exp(A::AbstractMatrix)
    F = eigen(A)
    V = F.vectors
    D = Diagonal(exp.(Complex.(F.values)))
    return real.(V * D * inv(V))
end

"""
    InfiniteHomographyHomotopy

A parameter homotopy that interpolates lines across views related by a homography H_∞.

For each line ℓ₀ passing through vanishing point v₀, and target line ℓ₁ passing through v₁ = H_∞ * v₀,
the homotopy interpolates using pencil coordinates:
1. Express ℓ₀ = λ₀*e₀ + μ₀*f₀ where (e₀, f₀) is a basis of lines through v₀
2. Transform basis: e_t = H_t^{-T} * e₀, f_t = H_t^{-T} * f₀ where H_t = exp(t*log(H_∞))
3. Interpolate: (λ_t, μ_t) = (1-t)*(λ₀, μ₀) + t*(λ₁, μ₁)
4. Reconstruct: ℓ_t = λ_t*e_t + μ_t*f_t

This ensures ℓ_t always passes through v_t = H_t * v₀.
"""
struct InfiniteHomographyHomotopy{T<:AbstractSystem} <: AbstractHomotopy
    F::T
    p::Vector{ComplexF64}  # start parameters (flattened lines)
    q::Vector{ComplexF64}  # target parameters (flattened lines)

    # Precomputed homography data
    H_inf::Matrix{Float64}        # H_∞
    H_inf_inv::Matrix{Float64}    # H_∞⁻¹
    H_inf_invT::Matrix{Float64}   # H_∞^{-T} = (H_∞⁻¹)ᵀ
    log_H_inf::Matrix{Float64}    # log(H_∞) for computing H_t = exp(t*log(H_∞))

    # Per-line precomputed data
    vanishing_points::Vector{Vector{Float64}}  # v₀ for each line
    pencil_bases_e::Vector{Vector{Float64}}    # e₀ for each line
    pencil_bases_f::Vector{Vector{Float64}}    # f₀ for each line
    lambda_start::Vector{Float64}              # λ₀ for each line
    mu_start::Vector{Float64}                  # μ₀ for each line
    lambda_target::Vector{Float64}             # λ₁ for each line
    mu_target::Vector{Float64}                 # μ₁ for each line

    # Cache
    t_cache::Base.RefValue{ComplexF64}
    pt::Vector{ComplexF64}
    taylor_pt::TaylorVector{2,ComplexF64}
end

function InfiniteHomographyHomotopy(
    F;
    start_parameters::AbstractVector,
    target_parameters::AbstractVector,
    H_inf::AbstractMatrix,
    vanishing_points::AbstractVector,
)
    InfiniteHomographyHomotopy(F, start_parameters, target_parameters, H_inf, vanishing_points)
end

function InfiniteHomographyHomotopy(
    F::ModelKit.System,
    p::AbstractVector,
    q::AbstractVector,
    H_inf::AbstractMatrix,
    vanishing_points::AbstractVector;
    compile::Union{Bool,Symbol} = true,
)
    InfiniteHomographyHomotopy(fixed(F; compile = compile), p, q, H_inf, vanishing_points)
end

function InfiniteHomographyHomotopy(
    F::AbstractSystem,
    p::AbstractVector,
    q::AbstractVector,
    H_inf::AbstractMatrix,
    vanishing_points::AbstractVector
)
    @assert length(p) == length(q) == nparameters(F)

    p̂ = Vector{ComplexF64}(p)
    q̂ = Vector{ComplexF64}(q)
    taylor_pt = TaylorVector{2}(ComplexF64, length(q))
    pt = copy(p̂)

    number_of_lines = length(p) ÷ 3
    @assert length(vanishing_points) == number_of_lines "Need one vanishing point per line"

    # Precompute homography matrices
    H_inf_f = Matrix{Float64}(H_inf)
    H_inf_inv = inv(H_inf_f)
    H_inf_invT = H_inf_inv'
    log_H_inf = matrix_log(H_inf_f)

    # Precompute per-line data
    vps = Vector{Vector{Float64}}(undef, number_of_lines)
    bases_e = Vector{Vector{Float64}}(undef, number_of_lines)
    bases_f = Vector{Vector{Float64}}(undef, number_of_lines)
    λ_start = zeros(Float64, number_of_lines)
    μ_start = zeros(Float64, number_of_lines)
    λ_target = zeros(Float64, number_of_lines)
    μ_target = zeros(Float64, number_of_lines)

    for i in 1:number_of_lines
        idx = (i-1)*3 + 1
        line_start = real.(p[idx:idx+2])
        line_target = real.(q[idx:idx+2])
        v0 = Vector{Float64}(vanishing_points[i])

        # Store vanishing point
        vps[i] = normalize(v0)

        # Compute pencil basis at v₀
        e0, f0 = pencil_basis(v0)
        bases_e[i] = e0
        bases_f[i] = f0

        # Express start line in pencil basis
        λ_start[i], μ_start[i] = express_in_pencil_basis(line_start, e0, f0)

        # Transform basis to target frame: e₁ = H_∞^{-T} * e₀
        e1 = normalize(H_inf_invT * e0)
        f1 = normalize(H_inf_invT * f0)

        # Express target line in transformed pencil basis
        λ_target[i], μ_target[i] = express_in_pencil_basis(line_target, e1, f1)

        # Normalize pencil coordinates to unit norm
        norm_start = sqrt(λ_start[i]^2 + μ_start[i]^2)
        λ_start[i] /= norm_start
        μ_start[i] /= norm_start

        norm_target = sqrt(λ_target[i]^2 + μ_target[i]^2)
        λ_target[i] /= norm_target
        μ_target[i] /= norm_target

        # Sign consistency: ensure we take the short path in projective space
        if λ_start[i] * λ_target[i] + μ_start[i] * μ_target[i] < 0
            λ_target[i] = -λ_target[i]
            μ_target[i] = -μ_target[i]
        end
    end

    InfiniteHomographyHomotopy(
        F,
        p̂,
        q̂,
        H_inf_f,
        H_inf_inv,
        H_inf_invT,
        log_H_inf,
        vps,
        bases_e,
        bases_f,
        λ_start,
        μ_start,
        λ_target,
        μ_target,
        Ref(complex(NaN)),
        pt,
        taylor_pt
    )
end

Base.size(H::InfiniteHomographyHomotopy) = size(H.F)

function start_parameters!(H::InfiniteHomographyHomotopy, p)
    H.p .= p
    H.t_cache[] = NaN
    H
end

function target_parameters!(H::InfiniteHomographyHomotopy, q)
    H.q .= q
    H.t_cache[] = NaN
    H
end

function parameters!(H::InfiniteHomographyHomotopy, p, q)
    H.p .= p
    H.q .= q
    H.t_cache[] = NaN
    H
end

"""
    interpolate_line(H, line_idx, t)

Interpolate line `line_idx` at parameter t ∈ [0,1].
Returns the interpolated line in homogeneous coordinates.
"""
function interpolate_line(H::InfiniteHomographyHomotopy, line_idx::Int, t::Real)
    # Compute H_t = exp(t * log(H_∞))
    H_t = matrix_exp(t * H.log_H_inf)
    H_t_invT = inv(H_t)'

    # Transform basis vectors: e_t = H_t^{-T} * e₀, f_t = H_t^{-T} * f₀
    e_t = H_t_invT * H.pencil_bases_e[line_idx]
    f_t = H_t_invT * H.pencil_bases_f[line_idx]

    # Interpolate pencil coordinates
    λ_t = (1 - t) * H.lambda_start[line_idx] + t * H.lambda_target[line_idx]
    μ_t = (1 - t) * H.mu_start[line_idx] + t * H.mu_target[line_idx]

    # Reconstruct interpolated line
    line_t = λ_t * e_t + μ_t * f_t

    return normalize(line_t)
end

"""
    tp!(H, t)

Compute interpolated parameters at time t and update the Taylor vector cache.
"""
function tp!(H::InfiniteHomographyHomotopy, tinput::Union{ComplexF64,Float64})
    tinput == H.t_cache[] && return H.taylor_pt
    t = real(tinput)

    number_of_lines = length(H.vanishing_points)
    parameters = zeros(3 * number_of_lines)

    # Compute H_t = exp(t * log(H_∞)) once for all lines
    H_t = matrix_exp(t * H.log_H_inf)
    H_t_invT = inv(H_t)'

    for i in 1:number_of_lines
        idx = (i-1)*3 + 1

        # Transform basis vectors: e_t = H_t^{-T} * e₀, f_t = H_t^{-T} * f₀
        e_t = H_t_invT * H.pencil_bases_e[i]
        f_t = H_t_invT * H.pencil_bases_f[i]

        # Interpolate pencil coordinates
        λ_t = (1 - t) * H.lambda_start[i] + t * H.lambda_target[i]
        μ_t = (1 - t) * H.mu_start[i] + t * H.mu_target[i]

        # Reconstruct interpolated line
        line_t = λ_t * e_t + μ_t * f_t

        # Normalize for numerical stability
        line_t = line_t / norm(line_t)

        parameters[idx:idx+2] = line_t
    end

    @inbounds for i = 1:length(H.taylor_pt)
        ptᵢ = parameters[i]
        H.pt[i] = ptᵢ
        H.taylor_pt[i] = (ptᵢ, H.p[i] - H.q[i])
    end
    H.t_cache[] = tinput

    H.taylor_pt
end

function ModelKit.evaluate!(u, H::InfiniteHomographyHomotopy, x, t)
    tp!(H, t)
    evaluate!(u, H.F, x, H.pt)
end

function ModelKit.evaluate_and_jacobian!(u, U, H::InfiniteHomographyHomotopy, x, t)
    tp!(H, t)
    evaluate_and_jacobian!(u, U, H.F, x, H.pt)
end

function ModelKit.taylor!(u, v::Val, H::InfiniteHomographyHomotopy, tx, t)
    taylor!(u, v, H.F, tx, tp!(H, t))
    u
end

"""
    verify_line_through_vanishing_point(H, line_idx, t)

Verify that the interpolated line at time t passes through the interpolated vanishing point.
Returns the absolute value of ℓ_tᵀ * v_t (should be ≈ 0).
"""
function verify_line_through_vanishing_point(H::InfiniteHomographyHomotopy, line_idx::Int, t::Real)
    # Compute v_t = H_t * v₀
    H_t = matrix_exp(t * H.log_H_inf)
    v_t = H_t * H.vanishing_points[line_idx]
    v_t = normalize(v_t)

    # Get interpolated line
    line_t = interpolate_line(H, line_idx, t)

    # Check incidence: ℓ_tᵀ * v_t should be 0
    return abs(dot(line_t, v_t))
end
