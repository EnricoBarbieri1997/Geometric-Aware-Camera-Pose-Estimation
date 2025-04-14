import json
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import os

# Label maps
method_labels = {
    "ours": "Ours",
    "quadric_based": "Quadric-based",
    "right_cylinder": "Right Cylinder",
    "zhang_4": "Zhang 4 views",
    "zhang_30": "Zhang 30 views"
}

metric_labels = {
    "delta_f": r"$\Delta f$",
    "delta_uv": r"$\Delta c$",
    "delta_skew": r"$\Delta$ $\gamma$",
}

method_supports = {
    "ours": {"delta_f": True, "delta_uv": True, "delta_skew": True},
    "quadric_based": {"delta_f": True, "delta_uv": True, "delta_skew": True},
    "right_cylinder": {"delta_f": True, "delta_uv": True, "delta_skew": False},
    "zhang_4": {"delta_f": True, "delta_uv": True, "delta_skew": False},
    "zhang_30": {"delta_f": True, "delta_uv": True, "delta_skew": False}
}

methods = list(method_labels.keys())
metrics = list(metric_labels.keys())
colors = ["blue", "purple", "green", "orange", "red"]
linestyles = ["-", "--", "-.", ":", ":"]

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

# Generate LaTeX
latex = []

latex.append(r"\begin{table}[H]")
latex.append(r"\centering")
latex.append(r"\begin{tabular}{l c l l l}")
latex.append(r"\toprule")
latex.append(r"Noise & Method & " +
             " & ".join(metric_labels[m] for m in metrics) + r" \\")
latex.append(r"\midrule")

for method in methods:
    noise_levels = sorted(grouped[metric][method].keys())
    for noise in noise_levels:
        supports = method_supports[method]
        row = [f"{noise:.2f}", method_labels[method]]

        for metric in metrics:
            if supports.get(metric, False):
                values = grouped[metric][method][noise]
                if values:
                    mean = np.mean(values)
                    row.append(f"{mean:.3f}")
                else:
                    row.append(r"--")
            else:
                row.append(r"--")

        latex.append(" & ".join(row) + r" \\")
    latex.append(r"\midrule")

latex.append(r"\bottomrule")
latex.append(r"\end{tabular}")
latex.append(r"\label{tab:SyntheticNoise}")
latex.append(r"\caption{Comparison of methods across different noise levels.}")
latex.append(r"\end{table}")

# Output to file or print
with open("table.tex", "w") as f:
    for line in latex:
        f.write(line + "\n")
