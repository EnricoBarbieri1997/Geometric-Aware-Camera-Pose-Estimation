import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from collections import defaultdict
import os

def customPlotFun(y, pos):
    if y == 0.0001:
        return "≤ 0.0001"
    else:
        return format(y)

# Label maps
method_labels = {
    "ours": "Ours",
    "quadric_based": "Quadric-based",
    "right_cylinder": "Right Cylinder",
    # "zhang_4": "Zhang 4 views",
    # "zhang_30": "Zhang 30 views"
}

metric_labels = {
    "delta_f": r"$\Delta f$",
    "delta_uv": r"$\Delta c$",
    "delta_skew": r"$\Delta$ $\gamma$",
    "delta_r": r"$\Delta R$",
    "delta_t": r"$\Delta t$",
    "success_rate": "Success Rate",
}

method_supports = {
    "ours": {"delta_f": True, "delta_uv": True, "delta_skew": True, "delta_r": True, "delta_t": True, "success_rate": True},
    "quadric_based": {"delta_f": True, "delta_uv": False, "delta_skew": False, "delta_r": True, "delta_t": True, "success_rate": True},
    "right_cylinder": {"delta_f": True, "delta_uv": True, "delta_skew": False, "delta_r": False, "delta_t": False, "success_rate": True},
    "zhang_4": {"delta_f": True, "delta_uv": True, "delta_skew": False, "delta_r": True, "delta_t": True, "success_rate": True},
    "zhang_30": {"delta_f": True, "delta_uv": True, "delta_skew": False, "delta_r": True, "delta_t": True, "success_rate": True},
}

methods = list(method_labels.keys())
metrics = list(metric_labels.keys())
colors = ["blue", "purple", "green", "orange", "red"]
linestyles = ["-", "--", "-.", ":", ":"]

# Load data
with open("./synthetic/results.json") as f:
    data = json.load(f)

# Group by metric → method → noise → values
grouped = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
for entry in data:
    noise = entry["noise"]
    method = entry["method"]
    for metric in metrics:
        values = entry.get(metric, [])
        if isinstance(values, list):
            for v in values:
                if v is not None:
                    grouped[metric][method][noise].append(v)
        else:
            if values is not None:
                grouped[metric][method][noise].append(values)

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
                vals = [max(v, 0.0001) for v in vals]
                mean = np.mean(vals)
                variance = np.var(vals)
                min_val = np.min(vals)
                max_val = np.max(vals)
                means.append(mean)
                plt.errorbar(
                    noise, mean, yerr=np.sqrt(variance), fmt='o', color=colors[idx]
                )
                plt.vlines(noise, min_val, max_val, color=colors[idx], alpha=0.5)
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
    plt.yscale("log")
    plt.gca().yaxis.set_major_formatter(customPlotFun)
    plt.title(f"{metric_labels[metric]} vs Noise")
    plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
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
