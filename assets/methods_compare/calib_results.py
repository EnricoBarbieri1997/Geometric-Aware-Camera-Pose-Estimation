# calib_results.py
import json
import os
import glob

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

def create_single_noise_result(method_code, noise_value, delta_f, delta_uv, delta_skew):
    """
    Creates a result entry for one noise level.
    """
    delta_f = [delta_f] if not isinstance(delta_f, list) else delta_f
    delta_uv = [delta_uv] if not isinstance(delta_uv, list) else delta_uv
    delta_skew = [delta_skew] if not isinstance(delta_skew, list) else delta_skew
    return {
        "noise": noise_value,
        "method": method_code,
        "delta_f": delta_f,
        "delta_uv": delta_uv,
        "delta_skew": delta_skew
    }

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
