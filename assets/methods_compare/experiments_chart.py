import json
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import os

# Label maps
method_labels = {
    "ours": "Ours",
    "quadric": "Quadric-based",
    "cylinder": "Right Cylinder",
    "zhang": "Zhang"
}

metric_labels = {
    "delta_f": r"$\Delta f$",
    "delta_uv": r"$\Delta uv$",
    "delta_skew": r"$\Delta$ $\Gamma$",
}

method_supports = {
    "ours": {"delta_f": True, "delta_uv": True, "delta_skew": True},
    "quadric": {"delta_f": True, "delta_uv": True, "delta_skew": True},
    "cylinder": {"delta_f": True, "delta_uv": True, "delta_skew": False},
    "zhang": {"delta_f": True, "delta_uv": True, "delta_skew": False}
}

methods = list(method_labels.keys())
metrics = list(metric_labels.keys())
colors = ["blue", "orange", "green", "red"]
linestyles = ["-", "--", "-.", ":"]

# Load data
with open("results.json") as f:
    data = json.load(f)

# Group by metric → method → noise → values
grouped = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
for entry in data:
    noise = entry["noise"]
    method = entry["method"]
    for metric in metrics:
        for v in entry.get(metric, []):
            if v is not None:
                grouped[metric][method][noise].append(v)

# Plot each metric
for metric in metrics:
    plt.figure(figsize=(5, 4))

    for idx, method in enumerate(methods):
        # Skip if the method does not support this metric
        if not method_supports.get(method, {}).get(metric, False):
            continue

        noise_levels = sorted(grouped[metric][method].keys())
        means = []
        for noise in noise_levels:
            vals = grouped[metric][method][noise]
            if vals:
                means.append(np.mean(vals))
            else:
                means.append(np.nan)

        plt.plot(
            noise_levels,
            means,
            label=method_labels[method],
            color=colors[idx],
            linestyle=linestyles[idx],
            marker='o'
        )

    plt.xlabel("Noise Level")
    plt.ylabel(metric_labels[metric])
    plt.title(f"{metric_labels[metric]} vs Noise")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    os.makedirs("plots", exist_ok=True)
    plt.savefig(f"plots/{metric}.pdf")
