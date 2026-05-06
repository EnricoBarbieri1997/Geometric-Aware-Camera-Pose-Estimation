# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Julia research project implementing cylinder-based camera resectioning (calibration) using homotopy continuation. The core idea: given 3D cylinders in a scene and their observed 2D silhouette contour lines in one or more camera views, recover camera intrinsic and extrinsic parameters by solving polynomial equation systems with `HomotopyContinuation.jl`.

The package is `CylindersBasedCameraResectioning` (Julia 1.11, `Project.toml` at root).

## Commands

### Running experiments (Julia REPL)

```julia
using CylindersBasedCameraResectioning

# Synthetic data experiment (main entry point)
CylindersBasedCameraResectioning.main()

# Real-scene experiments
CylindersBasedCameraResectioning.roller_coaster()
CylindersBasedCameraResectioning.hot_pipes()

# Benchmarking / lab comparisons
CylindersBasedCameraResectioning.Lab.compare_parameter_homotopies()
CylindersBasedCameraResectioning.Lab.compare_execution_times()

# Batch synthetic tests → ./tmp/reports/
CylindersBasedCameraResectioning.Report.multiple_seeds_multiple_configuration()
CylindersBasedCameraResectioning.Report.merge_reports("./tmp/reports/*", "./tmp/reports_merged.jls")
CylindersBasedCameraResectioning.Report.report_error_analysis("./tmp/reports_merged.jls", collect(0.0:0.0005:0.04))
```

### Tests

```bash
# Run all tests
julia --project=. -e 'using Pkg; Pkg.test()'

# Or from the REPL
]test
```

Individual test files are in `test/` and map to source modules (e.g. `test/geometry_tests.jl` → `src/geometry.jl`).

### Docker (headless / CI)

```bash
# Dev container (interactive Julia REPL, live-mounts src/)
docker compose -f docker-compose.dev.yml up

# Production run
docker compose up
```

The dev compose file mounts `src/`, `test/`, `assets/`, `tmp/` as volumes and drops into a Julia REPL. `Revise.jl` is auto-loaded via `~/.julia/config/startup.jl`.

### Nix dev shell

```bash
nix develop   # provides Julia 1.11, gcc, python3; sets LD_LIBRARY_PATH for GL libs
```

### Scripts

```bash
# Run batch noise experiments
./scripts/run.sh [seed_index] [noise_value] [output_path]

# Run noise sweep in parallel (PowerShell)
./scripts/run_noise.ps1
```

### Environment flags

| Variable | Default | Effect |
|---|---|---|
| `GUI_ENABLED` | `"true"` | Enables GLMakie visualisation; set `"false"` for headless runs |
| `ASSERTS_ENABLED` | `"false"` | Enables internal geometric assertion checks |

## Architecture

### Module hierarchy

All code lives under the `CylindersBasedCameraResectioning` module (entry: `src/CylindersBasedCameraResectioning.jl`). Sub-modules are included via `src/includes.jl` in dependency order:

```
Utils          – toleranced ≃ operator, random helpers, line format converters
Space          – 3-D primitive types (RotDeg), rotation / transformation factories
Camera         – CameraProperties, intrinsic matrix builder, random camera factory
Cylinder       – CylinderProperties, quadric/dual-quadric construction, calibration rigs
Geometry       – Line/Plane/Circle types, silhouette extraction (get_view), homogeneous math
Homotopies     – GeometricHomotopy and ParameterHomotopy wrappers around HomotopyContinuation
EquationSystems – Polynomial system builders; Problems sub-module holds
                  CylinderCameraContoursProblem and the IntrinsicParameters configuration enum
Scene          – SceneData / InstanceConfiguration, create_scene_instances_and_problems,
                 solver orchestration (intrinsic_rotation_system_setup, best_overall_solution_by_steps!)
Report         – Batch multi-seed experiment runner and error analysis
Lab            – Ad-hoc timing and comparison experiments
IO             – JSON / COLMAP / point-cloud readers (z-up convention)
Plotting       – GLMakie wrappers; GUI_ENABLED guards all display calls
Printing       – Camera error metrics, PrettyTables output
Conic, Formulas, Debug, Minimality, SolutionGenerator – supporting math modules
```

### Solving pipeline

1. **Scene construction** – `create_scene_instances_and_problems` generates (or loads) cylinders and per-view `CylinderCameraContoursProblem` instances containing tangent contour lines, dual quadrics, and points at infinity.
2. **Equation system** – `intrinsic_rotation_system_setup` builds a `HomotopyContinuation.System` whose variables are unknown intrinsic parameters (controlled by `IntrinsicParametersConfigurations` bitmask enum) and quaternion-parametrised rotations; observed line coordinates enter as parameters.
3. **Solve** – `solver_startsolutions` with `:total_degree` start system, or a `GeometricHomotopy` / `ParameterHomotopy` when start solutions from monodromy are available. Large solution counts are streamed in chunks of 500 000.
4. **Solution selection** – `best_overall_solution_by_steps!` filters real solutions, back-solves translation, and picks the minimum reprojection error camera.

### Intrinsic parameter configurations

`EquationSystems.Problems.IntrinsicParameters.Configurations` is a bitmask enum. Common values: `none`, `fₓ_fᵧ`, `fₓ_fᵧ_cₓ_cᵧ`, `fₓ_fᵧ_skew_cₓ_cᵧ`. The number of unknowns (and thus polynomial system size) grows with the configuration.

### Data formats

- **Scenes / views** are stored as JSON files under `assets/test_scenes/<scene>/`. `scene.json` contains 3-D cylinder geometry and ground-truth camera matrices; `views.json` contains observed 2-D contour lines per view.
- **Serialised results** use Julia's `Serialization` module (`.jls` / `.jld` files) in `tmp/`.
- **Line representation**: homogeneous 2-D lines stored as `Array{Float64, 3}` of shape `(n_cylinders, 2, 3)` — two tangent lines per cylinder.

### Assets
- `assets/test_scenes/` – synthetic and real scenes with cylinder geometry and observed contours. For each a `scene.json` (3-D geometry) and `views.json` (2-D contours) are stored.
- `assets/methods_compare/` - suite of Python scripts to run external methods for comparison, and create charts and table from data for visualization.
  - `zhang.py` - Zhang's method with OpenCV. Results in `synthetic/zhang_results.json`.
  - `right_cylinder.py` - right-cylinder method custom implementation. Results in `synthetic/right_cylinder_results.json`.
  - `experiments_charts_table.py` - reads results from JSON files and creates charts and tables for the paper. Outputs `synthetic/experiments_charts_table.png` and `synthetic/experiments_charts_table.tex`.
- `assets/article/` - A LaTeX article associated with the project. It should not be used for coding. It's mainly edited when content changes are needed.