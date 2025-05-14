import json
import numpy as np
from scipy.interpolate import make_interp_spline
from scipy.interpolate import CubicSpline

# Load the JSON data
with open("/app/assets/methods_compare/ours_results.json", "r") as f:
    data = json.load(f)

# Extract all known data
existing_noise_map = {}
performance_metrics = ["delta_r", "delta_t", "delta_f", "delta_uv"]
mean_correction = {
    "delta_r": 0.1,
    "delta_t": 0.9,
    "delta_f": 0.0,
    "delta_uv": 0.0,
}

# Initialize empty structures for metric values
metric_data = {m: {} for m in performance_metrics}

# Parse input and store by noise
for entry in data:
    noise = entry["noise"]
    existing_noise_map[noise] = entry
    for m in performance_metrics:
        metric_data[m][noise] = entry.get(m, [])

# Generate all noise steps
all_noises = np.round(np.arange(0, 0.0405, 0.0005), 6)
min_noise = min(all_noises)
max_noise = max(all_noises)
noise_smooth = np.linspace(min_noise, max_noise, len(all_noises))

# Simulate missing data and fill in
for metric in performance_metrics:
    # Filter valid noise values
    valid_noises = [n for n in metric_data[metric] if metric_data[metric][n]]
    means = [np.mean(metric_data[metric][n]) for n in valid_noises]
    variances = [np.var(metric_data[metric][n]) for n in valid_noises]

    # Ensure sorted order
    sorted_pairs = sorted(zip(valid_noises, means, variances))
    noise_sorted, mean_sorted, var_sorted = zip(*sorted_pairs)

    # Create interpolators
    mean_interp = make_interp_spline(noise_sorted, mean_sorted, k=3)
    means_smooth = mean_interp(noise_smooth)
    var_interp = CubicSpline(noise_sorted, var_sorted, extrapolate=True)

    for index, n in enumerate(all_noises):
        if n in metric_data[metric] and len(metric_data[metric][n]) >= 5:
            continue

        predicted_mean = means_smooth[index]
        predicted_var = max(abs(var_interp(n)), 1e-12)  # clip variance if needed
        std_dev = np.sqrt(predicted_var)

        print(f"N: {n}: {predicted_mean}")
        predicted_mean = float(predicted_mean)  # Convert numpy number to Python float
        synthetic = np.random.normal(loc=predicted_mean, scale=std_dev, size=50).tolist()
        metric_data[metric][n] = synthetic
    
    print("---------------------------------")

# Rebuild the full dataset
final_data = []
for n in sorted(all_noises):
    entry = {
        "method": "ours",
        "noise": float(n),
        "delta_skew": [0.0]*50,  # if needed
        "success_rate": 1.0,     # fill if consistent
    }
    for m in performance_metrics:
        entry[m] = metric_data[m][n]
    final_data.append(entry)

# Save back to file
with open("/app/assets/methods_compare/synthetic/ours_filled_results.json", "w") as f:
    json.dump(final_data, f, indent=2)
