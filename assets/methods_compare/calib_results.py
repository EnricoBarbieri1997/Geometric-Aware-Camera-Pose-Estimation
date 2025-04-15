# calib_results.py
import json
import os
import glob
import numpy as np

def iterations_results_to_metrics(iterations_results):
    delta_f = [res[0] for res in iterations_results]
    delta_uv = [res[1] for res in iterations_results]
    delta_skew = [res[2] for res in iterations_results]
    success_rate = sum(res[3] for res in iterations_results)
    delta_r = None
    delta_t = None
    if len(iterations_results[0]) > 4:
        delta_r = [res[4] for res in iterations_results]
    if len(iterations_results[0]) > 5:
        delta_t = [res[5] for res in iterations_results]

    return delta_f, delta_uv, delta_skew, success_rate, delta_r, delta_t

def generate_K():
    # Intrinsic matrix with noise
    fx = 3500 + np.random.normal(0, 50)
    fy = 3500 + np.random.normal(0, 50)
    cx = 1750 + np.random.normal(0, 10)
    cy = 1750 + np.random.normal(0, 10)
    K = np.array([
        [fx,  0, cx],
        [0,  fy, cy],
        [0,   0,  1]
    ])
    return K

def normalized_diff(calculated, truth):
    if calculated == 0 and truth == 0:
        return 0.0
    denominator = truth if truth != 0 else calculated
    return abs(calculated - truth) / abs(denominator)

def intrinsic_difference(calculated, truth):
    truth = truth / truth[2, 2]
    calculated = calculated / calculated[2, 2]
    fx, fy = calculated[0, 0], calculated[1, 1]
    cx, cy = calculated[0, 2], calculated[1, 2]
    skew = calculated[0, 1]

    fx_t, fy_t = truth[0, 0], truth[1, 1]
    cx_t, cy_t = truth[0, 2], truth[1, 2]
    skew_t = truth[0, 1]

    deltaF = (normalized_diff(fx, fx_t) + normalized_diff(fy, fy_t)) / 2
    deltaUV = (normalized_diff(cx, cx_t) + normalized_diff(cy, cy_t)) / 2
    deltaSkew = normalized_diff(skew, skew_t)

    return [deltaF, deltaUV, deltaSkew]

def rotations_difference(R1: np.ndarray, R2: np.ndarray) -> float:
    trace_val = np.trace(R1 @ R2.T)
    return np.degrees(np.arccos((trace_val - 1) / 2))

def translations_difference(t1: np.ndarray, t2: np.ndarray) -> float:
    return np.linalg.norm(t1 - t2)

def create_single_noise_result(method_code, noise_value, delta_f, delta_uv, delta_skew, success_rate, delta_r=None, delta_t=None):
    """
    Creates a result entry for one noise level.
    """
    delta_f = [delta_f] if not isinstance(delta_f, list) else delta_f
    delta_uv = [delta_uv] if not isinstance(delta_uv, list) else delta_uv
    delta_skew = [delta_skew] if not isinstance(delta_skew, list) else delta_skew

    metrics = {
        "noise": noise_value,
        "method": method_code,
        "delta_f": delta_f,
        "delta_uv": delta_uv,
        "delta_skew": delta_skew,
        "success_rate": success_rate,
    }
    if delta_r is not None:
        delta_r = [delta_r] if not isinstance(delta_r, list) else delta_r
        metrics["delta_r"] = delta_r
    if delta_t is not None:
        delta_t = [delta_t] if not isinstance(delta_t, list) else delta_t
        metrics["delta_t"] = delta_t

    return metrics

def save_results_to_json(filename, results):
    """
    Saves a list of result entries to a JSON file.
    """
    with open(filename, "w") as f:
        json.dump(results, f, indent=2)

def merge_all_results_json(output_file="results.json", pattern="*_results.json"):
    """
    Merges all *_results.json files in the current folder into a single file.
    """
    merged = []
    for file in glob.glob(pattern):
        with open(file) as f:
            try:
                merged.extend(json.load(f))
            except Exception as e:
                print(f"⚠️ Skipping {file} due to error: {e}")
    with open(output_file, "w") as f:
        json.dump(merged, f, indent=2)
    print(f"✅ Merged {len(merged)} entries into {output_file}")

# Allow running this script directly to merge
if __name__ == "__main__":
    merge_all_results_json()
