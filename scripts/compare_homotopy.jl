#!/usr/bin/env julia
"""
Compare geometric-parametric homotopy against the Ceres solver results.

Loads the JSON scene file produced by generate_ceres_scene.jl, reconstructs
CylinderCameraContoursProblem objects, runs monodromy + GeometricHomotopy and
prints the same table format as the Ceres solver.  Falls back to the standard
straight-line ParameterHomotopy if the geometric path fails.

The polynomial system solved is the gradient of Σ(l'KRd)² — the same
Gauss-Helmert algebraic cost that Ceres minimises — set to zero.  Solving its
critical-point system with homotopy continuation finds the global minimum and
all other critical points simultaneously.  Because a minimum always exists in ℝ,
this formulation works with both noise-free and noisy silhouette lines.

Usage:
  julia --project=. scripts/compare_homotopy.jl
  julia --project=. scripts/compare_homotopy.jl ceres_solver/scene_data.json
  julia --project=. scripts/compare_homotopy.jl ceres_solver/scene_data_noisy.json
"""

ENV["GUI_ENABLED"] = "false"
using Pkg
Pkg.activate(dirname(dirname(abspath(@__FILE__))))

using CylindersBasedCameraResectioning
using HomotopyContinuation
using LinearAlgebra
using JSON
using Rotations

const CamM   = CylindersBasedCameraResectioning.Camera
const SceneM = CylindersBasedCameraResectioning.Scene
const HomotM = CylindersBasedCameraResectioning.Homotopies

const EquSysProblems = CylindersBasedCameraResectioning.EquationSystems.Problems
const IntrinsicConf  = EquSysProblems.IntrinsicParameters.Configurations

# ─── CLI arguments ─────────────────────────────────────────────────────────────

json_path = length(ARGS) > 0 ? ARGS[1] : "ceres_solver/scene_data.json"
println("Loading scene from $json_path")
data = JSON.parsefile(json_path)

n_views     = data["n_views"]
n_cylinders = data["n_cylinders"]

intrinsic_configuration = IntrinsicConf.fₓ_fᵧ_cₓ_cᵧ

# ─── Group residuals by view ───────────────────────────────────────────────────

residuals_by_view = [Any[] for _ in 1:n_views]
for entry in data["residuals"]
    v = entry["view"] + 1   # 1-indexed
    push!(residuals_by_view[v], entry)
end

# ─── Build CylinderCameraContoursProblem from JSON entries ────────────────────

function orient_line(l)
    # Match the sign convention used by problems_from_scene: divide by the
    # third component first so that c is always rendered positive, then L2-normalize.
    # This is critical for GeometricHomotopy: if start and target lines carry
    # opposite signs the "angle between" them is ≈π instead of ≈0 and the
    # homotopy path diverges.
    l = l[3] < 0 ? -l : l
    return normalize(l)
end

function build_problem(entries, use_start_lines)
    n = length(entries)
    lines    = Matrix{Float64}(undef, n, 3)
    pts_inf  = Matrix{Float64}(undef, n, 3)
    for (i, e) in enumerate(entries)
        raw_l = use_start_lines ? Float64.(e["start_line"]) : Float64.(e["line"])
        lines[i, :]   = orient_line(raw_l)
        pts_inf[i, :] = normalize(Float64.(e["axis"]))
    end

    cam = CamM.CameraProperties()
    cam.intrinsic = [2600.0 0.0 960.0; 0.0 1500.0 540.0; 0.0 0.0 1.0]

    validation = EquSysProblems.CylinderCameraContoursProblemValidationData(
        Matrix{Float64}(undef, 0, 3),
        Matrix{Float64}(undef, 0, 3),
        Array{Float64}(undef, 0, 4, 4),
    )

    return EquSysProblems.CylinderCameraContoursProblem(
        cam, lines, lines, pts_inf,
        Array{Float64}(undef, n, 4, 4),
        validation,
        UInt8(intrinsic_configuration),
    )
end

target_problems = [build_problem(residuals_by_view[v], false) for v in 1:n_views]
start_problems  = [build_problem(residuals_by_view[v], true)  for v in 1:n_views]

# ─── Build polynomial system and parameters ────────────────────────────────────

rotation_intrinsic_system, start_params = SceneM.intrinsic_rotation_system_setup(
    start_problems; minimization = true, intrinsic_configuration,
)
_, target_params = SceneM.intrinsic_rotation_system_setup(
    target_problems; minimization = true, intrinsic_configuration,
)

# ─── Helpers: rotation matrix ↔ Cayley / angle-axis ───────────────────────────

function aa_to_rotmat(aa_vec)
    θ = norm(aa_vec)
    θ < 1e-10 && return Matrix{Float64}(I, 3, 3)
    return Matrix{Float64}(RotMatrix(AngleAxis(θ, (aa_vec / θ)...)))
end

function rotmat_to_cayley(R_mat)
    quat = Rotations.params(QuatRotation(RotMatrix{3}(R_mat)))
    return quat[2:4] / quat[1]
end

function build_rotation_matrix_normalized(x, y, z)
    k = 1 + x^2 + y^2 + z^2
    return (1/k) * [
        1+x^2-y^2-z^2  2*x*y-2*z    2*y+2*x*z;
        2*z+2*x*y      1-x^2+y^2-z^2  2*y*z-2*x;
        2*x*z-2*y      2*x+2*y*z    1-x^2-y^2+z^2
    ]
end

# ─── Known start solution: perturbed intrinsics + Cayley params ───────────────

factor = 1.0 / 3000.0

sintr = data["start"]["intrinsics"]
start_fx = Float64(sintr["fx"])
start_fy = Float64(sintr["fy"])
start_cx = Float64(sintr["cx"])
start_cy = Float64(sintr["cy"])

start_cayley = [
    rotmat_to_cayley(aa_to_rotmat(Float64.(data["start"]["rotations"][v]["angle_axis"])))
    for v in 1:n_views
]

known_solution = [
    start_fx * factor, start_fy * factor,
    start_cx * factor, start_cy * factor,
    vcat(start_cayley...)...,
]

# ─── Geometric homotopy: track known start solution → target ──────────────────
#
# For the minimization system (∇Σ(l'KRd)² = 0) a global minimum always exists
# in ℝ, so we can track the single known start solution directly — no monodromy
# required.  The geometric homotopy rotates each line from its start direction to
# its target direction, keeping the path smooth and real throughout.

println("\nTracking known solution via GeometricHomotopy (minimization system) …")
geo_result = solve(
    HomotM.GeometricHomotopy(
        rotation_intrinsic_system,
        start_parameters = start_params,
        target_parameters = target_params,
        compile = false,   # skip LLVM compilation; degree-5 system is too large to JIT
    ),
    [known_solution];
    show_progress = true,
)
println("GeometricHomotopy: $(nsolutions(geo_result)) non-singular, $(nreal(geo_result)) real, $(nfailed(geo_result)) failed")

# ─── Straight-line parameter homotopy fallback ────────────────────────────────

lin_result = nothing
if nreal(geo_result) == 0
    println("\nGeometric path failed — running straight-line ParameterHomotopy …")
    lin_result = solve(
        rotation_intrinsic_system,
        [known_solution];
        start_parameters = start_params,
        target_parameters = target_params,
        compile = false,
        show_progress = true,
    )
    println("ParameterHomotopy: $(nsolutions(lin_result)) non-singular, $(nreal(lin_result)) real, $(nfailed(lin_result)) failed")
end

# Accept solutions with small imaginary/real ratio
function approx_real_solutions(result; tol=1e-4)
    sols = solutions(result)
    isempty(sols) && return Vector{Float64}[]
    return [real.(s) for s in sols if norm(imag.(s)) < tol * (norm(real.(s)) + 1.0)]
end

geo_approx = approx_real_solutions(geo_result; tol=1e-4)
lin_approx = !isnothing(lin_result) ? approx_real_solutions(lin_result; tol=1e-4) : Vector{Float64}[]

println("\nApprox-real solutions (imag/real tol=1e-4): geo=$(length(geo_approx)), lin=$(length(lin_approx))")

# Choose the best result
use_sols, method_name = if !isempty(geo_approx)
    geo_approx, "GeometricHomotopy (minimization)"
elseif !isempty(lin_approx)
    lin_approx, "ParameterHomotopy (minimization, fallback)"
elseif nreal(geo_result) > 0
    real_solutions(geo_result), "GeometricHomotopy (minimization)"
elseif !isnothing(lin_result) && nreal(lin_result) > 0
    real_solutions(lin_result), "ParameterHomotopy (minimization, fallback)"
else
    println("\nNo real solutions from either method — cannot build table.")
    exit(1)
end

println("\nUsing results from: $method_name")

# ─── True solution in scaled Cayley space ─────────────────────────────────────

tintr = data["true"]["intrinsics"]
true_fx = Float64(tintr["fx"])
true_fy = Float64(tintr["fy"])
true_cx = Float64(tintr["cx"])
true_cy = Float64(tintr["cy"])

true_cayley = [
    rotmat_to_cayley(aa_to_rotmat(Float64.(data["true"]["rotations"][v]["angle_axis"])))
    for v in 1:n_views
]

true_solution_scaled = [
    true_fx * factor, true_fy * factor,
    true_cx * factor, true_cy * factor,
    vcat(true_cayley...)...,
]

best_idx = findmin([norm(true_solution_scaled - s) for s in use_sols])[2]
best     = use_sols[best_idx]

solved_fx = best[1] / factor
solved_fy = best[2] / factor
solved_cx = best[3] / factor
solved_cy = best[4] / factor

# ─── Print table ──────────────────────────────────────────────────────────────

rel_err(solved, truth) = truth == 0.0 ? abs(solved) : abs((solved - truth) / truth) * 100.0

W = 10
fmt(x) = lpad(round(x, digits=4), W)

println("\n" * "═"^59)
println("  Intrinsic parameters  ($method_name)")
println("─"^59)
println(lpad("param",W), lpad("true",W), lpad("start",W), lpad("solved",W), lpad("err (%)",W))
println("─"^59)
for (name, truth, start, solved) in [
    ("fx", true_fx, start_fx, solved_fx),
    ("fy", true_fy, start_fy, solved_fy),
    ("cx", true_cx, start_cx, solved_cx),
    ("cy", true_cy, start_cy, solved_cy),
]
    println(lpad(name,W), fmt(truth), fmt(start), fmt(solved),
            lpad(string(round(rel_err(solved,truth), digits=4))*"%", W))
end

function rotation_error_deg(R_true, R_solved)
    tr = sum(R_true .* R_solved)
    acosd(clamp((tr - 1.0) / 2.0, -1.0, 1.0))
end

for v in 1:n_views
    aa_true  = Float64.(data["true"]["rotations"][v]["angle_axis"])
    aa_start = Float64.(data["start"]["rotations"][v]["angle_axis"])

    R_true  = aa_to_rotmat(aa_true)
    R_start = aa_to_rotmat(aa_start)

    offset = 4 + (v - 1) * 3
    c1, c2, c3 = best[offset+1], best[offset+2], best[offset+3]
    R_solved = build_rotation_matrix_normalized(c1, c2, c3)

    tc = true_cayley[v]
    sc = start_cayley[v]

    println("─"^59)
    println("  Rotation (Cayley params)  -- view $(v-1)")
    println("─"^59)
    println(lpad("param",W), lpad("true",W), lpad("start",W), lpad("solved",W), lpad("|err|",W))
    println("─"^59)
    for (k_name, t, s, sol) in zip(["c1","c2","c3"], tc, sc, [c1,c2,c3])
        println(lpad(k_name,W), fmt(t), fmt(s), fmt(sol), fmt(abs(sol-t)))
    end

    ang_start  = rotation_error_deg(R_true, R_start)
    ang_solved = rotation_error_deg(R_true, R_solved)
    println("  Angular error:  start $(round(ang_start, digits=4)) deg   solved $(round(ang_solved, digits=4)) deg")
end
println("═"^59)
