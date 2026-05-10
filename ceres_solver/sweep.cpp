/**
 * Batch sweep: cylinder-based camera calibration via Ceres over
 * 81 noise levels × 50 trials.
 *
 * For each σ ∈ {0.0, 0.0005, …, 0.04} and 50 independent trials,
 * pixel noise is added to the observed line coordinates (matching
 * Julia's add_noise_to_line), Ceres is run from the same perturbed
 * starting point, and calibration errors are recorded.
 *
 * Input:  scene_data.json  (noiseless, from scripts/generate_ceres_scene.jl)
 * Output: ceres_results.json  written to the path given as argv[2]
 *         (default: ../../assets/methods_compare/synthetic/ceres_results.json)
 *
 * Build:  cmake .. && make calibrate_sweep
 * Run:    ./calibrate_sweep [scene_data.json] [output_path]
 */

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <nlohmann/json.hpp>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using json = nlohmann::json;
namespace fs = std::filesystem;

// ─── Cost functor (identical to main.cpp) ────────────────────────────────────

struct CylinderLineCost {
    CylinderLineCost(std::array<double,3> line, std::array<double,3> axis)
        : l_(line), d_(axis) {}

    template <typename T>
    bool operator()(const T* const intr, const T* const aa, T* residual) const {
        const T d_w[3] = {T(d_[0]), T(d_[1]), T(d_[2])};
        T d_c[3];
        ceres::AngleAxisRotatePoint(aa, d_w, d_c);
        T KRd[3];
        KRd[0] = intr[0] * d_c[0] + intr[2] * d_c[2];
        KRd[1] = intr[1] * d_c[1] + intr[3] * d_c[2];
        KRd[2] = d_c[2];
        residual[0] = T(l_[0])*KRd[0] + T(l_[1])*KRd[1] + T(l_[2])*KRd[2];
        return true;
    }

    static ceres::CostFunction* Create(std::array<double,3> line,
                                       std::array<double,3> axis) {
        return new ceres::AutoDiffCostFunction<CylinderLineCost, 1, 4, 3>(
            new CylinderLineCost(line, axis));
    }

    std::array<double,3> l_, d_;
};

// ─── Noise model ─────────────────────────────────────────────────────────────
//
// Mirrors Julia's Scene.add_noise_to_line:
//   p1 = [0,   intercept(0),   1] + normalize([rand-0.5, rand-0.5, 0]) * noise
//   p2 = [1000, intercept(1000), 1] + normalize([rand-0.5, rand-0.5, 0]) * noise
//   return normalize(cross(p1, p2))

static std::array<double,3> add_noise_to_line(
    const std::array<double,3>& l,
    double noise,
    std::mt19937& rng)
{
    std::uniform_real_distribution<double> ud(0.0, 1.0);

    auto noisy_pt = [&](double x) -> std::array<double,3> {
        double y = -(l[0]*x + l[2]) / l[1];
        double dx = ud(rng) - 0.5;
        double dy = ud(rng) - 0.5;
        double len = std::sqrt(dx*dx + dy*dy);
        if (len > 1e-12) { dx /= len; dy /= len; }
        return {x + dx*noise, y + dy*noise, 1.0};
    };

    auto p1 = noisy_pt(0.0);
    auto p2 = noisy_pt(1000.0);

    std::array<double,3> out = {
        p1[1]*p2[2] - p1[2]*p2[1],
        p1[2]*p2[0] - p1[0]*p2[2],
        p1[0]*p2[1] - p1[1]*p2[0]
    };
    double n = std::sqrt(out[0]*out[0] + out[1]*out[1] + out[2]*out[2]);
    if (n > 1e-12) { out[0] /= n; out[1] /= n; out[2] /= n; }
    return out;
}

// ─── Error metrics ────────────────────────────────────────────────────────────
//
// Matches Julia Utils: normalized_diff, intrinsic_difference, rotations_difference.

static double norm_diff(double solved, double truth) {
    if (truth == 0.0 && solved == 0.0) return 0.0;
    double denom = (truth != 0.0) ? std::abs(truth) : std::abs(solved);
    return std::abs(solved - truth) / denom;
}

static double rotation_error_deg(const double aa_true[3],
                                  const double aa_solved[3]) {
    double R_t[9], R_s[9];
    ceres::AngleAxisToRotationMatrix(aa_true,   R_t);
    ceres::AngleAxisToRotationMatrix(aa_solved, R_s);
    double tr = 0.0;
    for (int k = 0; k < 9; ++k) tr += R_t[k] * R_s[k];
    double cos_ang = std::max(-1.0, std::min(1.0, (tr - 1.0) / 2.0));
    return std::acos(cos_ang) * 180.0 / M_PI;
}

// ─── Main ─────────────────────────────────────────────────────────────────────

int main(int argc, char** argv) {
    // ── Argument parsing ──────────────────────────────────────────────────────
    // Usage: calibrate_sweep [scene.json] [out.json] [--random-start]
    std::string scene_path = "scene_data.json";
    std::string out_path   = "../../assets/methods_compare/synthetic/ceres_results.json";
    bool random_start = false;

    int positional = 0;
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg == "--random-start") {
            random_start = true;
        } else if (positional == 0) {
            scene_path = arg; ++positional;
        } else if (positional == 1) {
            out_path = arg; ++positional;
        }
    }

    std::cout << "start mode: " << (random_start ? "random" : "perturbed") << "\n";

    // ── Load base scene (noiseless) ───────────────────────────────────────────
    std::ifstream f(scene_path);
    if (!f.is_open()) {
        std::cerr << "Cannot open " << scene_path << "\n";
        return 1;
    }
    json data = json::parse(f);

    const int n_views = data["n_views"];

    const double true_fx = data["true"]["intrinsics"]["fx"];
    const double true_fy = data["true"]["intrinsics"]["fy"];
    const double true_cx = data["true"]["intrinsics"]["cx"];
    const double true_cy = data["true"]["intrinsics"]["cy"];

    std::vector<std::array<double,3>> aa_true(n_views);
    for (int v = 0; v < n_views; ++v) {
        auto& s = data["true"]["rotations"][v]["angle_axis"];
        aa_true[v] = {s[0], s[1], s[2]};
    }

    const double start_fx = data["start"]["intrinsics"]["fx"];
    const double start_fy = data["start"]["intrinsics"]["fy"];
    const double start_cx = data["start"]["intrinsics"]["cx"];
    const double start_cy = data["start"]["intrinsics"]["cy"];

    std::vector<std::array<double,3>> aa_start(n_views);
    for (int v = 0; v < n_views; ++v) {
        auto& s = data["start"]["rotations"][v]["angle_axis"];
        aa_start[v] = {s[0], s[1], s[2]};
    }

    struct ResidualEntry { int view; std::array<double,3> line, axis; };
    std::vector<ResidualEntry> base_residuals;
    for (auto& e : data["residuals"]) {
        base_residuals.push_back({
            e["view"],
            {e["line"][0], e["line"][1], e["line"][2]},
            {e["axis"][0], e["axis"][1], e["axis"][2]}
        });
    }

    // ── Sweep configuration ───────────────────────────────────────────────────
    constexpr int    N_TRIALS   = 50;
    constexpr double NOISE_MIN  = 0.0;
    constexpr double NOISE_MAX  = 0.04;
    constexpr double NOISE_STEP = 0.0005;
    const int N_NOISE =
        static_cast<int>(std::round((NOISE_MAX - NOISE_MIN) / NOISE_STEP)) + 1;

    // ── Solver options (silent) ───────────────────────────────────────────────
    ceres::Solver::Options opts;
    opts.linear_solver_type            = ceres::DENSE_QR;
    opts.minimizer_progress_to_stdout  = false;
    opts.max_num_iterations            = 500;
    opts.function_tolerance            = 1e-14;
    opts.gradient_tolerance            = 1e-14;
    opts.parameter_tolerance           = 1e-14;
    opts.logging_type                  = ceres::SILENT;

    // ── Run sweep ─────────────────────────────────────────────────────────────
    json results = json::array();

    for (int ni = 0; ni < N_NOISE; ++ni) {
        double noise = std::round((NOISE_MIN + ni * NOISE_STEP) * 100000.0) / 100000.0;

        std::vector<double> v_r, v_t, v_f, v_uv, v_skew;
        int n_success = 0;

        std::cout << "noise=" << std::fixed << std::setprecision(4) << noise
                  << "  [" << ni+1 << "/" << N_NOISE << "]  ";
        std::cout.flush();

        for (int trial = 0; trial < N_TRIALS; ++trial) {
            // Reproducible noise per trial (seed = trial+1)
            std::mt19937 rng(static_cast<uint32_t>(trial + 1));

            // Apply noise to line observations
            std::vector<std::array<double,3>> noisy_lines;
            noisy_lines.reserve(base_residuals.size());
            for (auto& re : base_residuals) {
                if (noise > 0.0)
                    noisy_lines.push_back(add_noise_to_line(re.line, noise, rng));
                else
                    noisy_lines.push_back(re.line);
            }

            // Starting parameters: either the stored perturbed solution or random
            double intr[4];
            std::vector<std::array<double,3>> aa(n_views);

            if (!random_start) {
                intr[0] = start_fx; intr[1] = start_fy;
                intr[2] = start_cx; intr[3] = start_cy;
                for (int v = 0; v < n_views; ++v) aa[v] = aa_start[v];
            } else {
                // Independent RNG so noise and start don't correlate
                std::mt19937 srng(static_cast<uint32_t>(trial + 10000));
                std::uniform_real_distribution<double> f_dist(500.0, 5000.0);
                std::uniform_real_distribution<double> cx_dist(0.0, 1920.0);
                std::uniform_real_distribution<double> cy_dist(0.0, 1080.0);
                intr[0] = f_dist(srng);
                intr[1] = f_dist(srng);
                intr[2] = cx_dist(srng);
                intr[3] = cy_dist(srng);

                std::normal_distribution<double> gauss(0.0, 1.0);
                std::uniform_real_distribution<double> angle_dist(0.0, M_PI);
                for (int v = 0; v < n_views; ++v) {
                    double ax = gauss(srng), ay = gauss(srng), az = gauss(srng);
                    double len = std::sqrt(ax*ax + ay*ay + az*az);
                    if (len > 1e-12) { ax /= len; ay /= len; az /= len; }
                    double angle = angle_dist(srng);
                    aa[v] = {ax*angle, ay*angle, az*angle};
                }
            }

            // Build and solve
            ceres::Problem problem;
            for (size_t ri = 0; ri < base_residuals.size(); ++ri) {
                int v = base_residuals[ri].view;
                problem.AddResidualBlock(
                    CylinderLineCost::Create(noisy_lines[ri], base_residuals[ri].axis),
                    nullptr,
                    intr, aa[v].data()
                );
            }
            ceres::Solver::Summary summary;
            ceres::Solve(opts, &problem, &summary);

            bool converged =
                summary.termination_type == ceres::CONVERGENCE ||
                summary.termination_type == ceres::USER_SUCCESS;

            // Metrics (match Julia Utils conventions)
            double delta_r = 0.0;
            for (int v = 0; v < n_views; ++v)
                delta_r += rotation_error_deg(aa_true[v].data(), aa[v].data());
            delta_r /= n_views;

            // delta_f: mean relative focal-length error (matches intrinsic_difference)
            double delta_f  = (norm_diff(intr[0], true_fx) +
                               norm_diff(intr[1], true_fy)) / 2.0;

            // delta_uv: mean relative principal-point error
            double delta_uv = (norm_diff(intr[2], true_cx) +
                               norm_diff(intr[3], true_cy)) / 2.0;

            if (converged) ++n_success;

            v_r.push_back(delta_r);
            v_t.push_back(0.0);
            v_f.push_back(delta_f);
            v_uv.push_back(delta_uv);
            v_skew.push_back(0.0);
        }

        double success_rate = static_cast<double>(n_success) / N_TRIALS;
        std::cout << "success=" << n_success << "/" << N_TRIALS << "\n";

        json entry;
        entry["noise"]        = noise;
        entry["method"]       = "ceres";
        entry["delta_r"]      = v_r;
        entry["delta_t"]      = v_t;
        entry["delta_f"]      = v_f;
        entry["delta_uv"]     = v_uv;
        entry["delta_skew"]   = v_skew;
        entry["success_rate"] = success_rate;
        results.push_back(std::move(entry));
    }

    // ── Write output ──────────────────────────────────────────────────────────
    fs::path out(out_path);
    if (out.has_parent_path())
        fs::create_directories(out.parent_path());

    std::ofstream ofs(out_path);
    if (!ofs.is_open()) {
        std::cerr << "Cannot write to " << out_path << "\n";
        return 1;
    }
    ofs << results.dump(2) << "\n";
    std::cout << "\nResults written to " << out_path << "\n";
    return 0;
}
