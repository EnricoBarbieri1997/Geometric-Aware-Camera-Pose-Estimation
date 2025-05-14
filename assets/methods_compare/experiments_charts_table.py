import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import LogLocator, FuncFormatter
from collections import defaultdict
import os
from scipy.interpolate import make_interp_spline

def customPlotFun(y, pos):
    if y <= 0.0001:
        return "≤ 0.0001"
    else:
        return f"{y:.5f}".rstrip('0').rstrip('.')

# Label maps
method_labels = {
    "ours": "Ours",
    "ours_skew": "Ours",
    "quadric_based": "Gummeson",
    "right_cylinder": "Ding",
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
    "ours": {"delta_f": True, "delta_uv": True, "delta_skew": False, "delta_r": True, "delta_t": True, "success_rate": True},
    "ours_skew": {"delta_f": False, "delta_uv": False, "delta_skew": True, "delta_r": False, "delta_t": False, "success_rate": False},
    "quadric_based": {"delta_f": False, "delta_uv": False, "delta_skew": False, "delta_r": True, "delta_t": True, "success_rate": True},
    "right_cylinder": {"delta_f": True, "delta_uv": True, "delta_skew": False, "delta_r": False, "delta_t": False, "success_rate": True},
    "zhang_4": {"delta_f": True, "delta_uv": True, "delta_skew": False, "delta_r": True, "delta_t": True, "success_rate": True},
    "zhang_30": {"delta_f": True, "delta_uv": True, "delta_skew": False, "delta_r": True, "delta_t": True, "success_rate": True},
}

methods = list(method_labels.keys())
metrics = list(metric_labels.keys())
colors = ["blue", "blue", "purple", "green", "orange", "red"]
linestyles = ["-", "-", "--", "-.", ":", ":"]

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
                vals = np.array([max(v, 0.0001) for v in vals])
                variance = np.var(vals)
                q25 = np.percentile(vals, 25)
                q75 = np.percentile(vals, 75)
                vals_iqr = vals[(vals >= q25) & (vals <= q75)]
                mean = np.mean(vals_iqr) if len(vals_iqr) > 0 else np.mean(vals)
                means.append(mean)
                err_low = max(mean - q25, 0)
                err_high = max(q75 - mean,0)

                if (err_low != 0 and err_high != 0):
                    plt.errorbar(
                        noise * 100, mean, yerr=[[err_low], [err_high]], fmt='o', color=colors[idx]
                    )
            else:
                means.append(np.nan)

        x_values = np.array(noise_levels) * 100
        if len(x_values) >= 4:
            x_smooth = np.linspace(x_values.min(), x_values.max(), 300)
            spline = make_interp_spline(x_values, means, k=3)  # k=3 → cubic
            y_smooth = spline(x_smooth)
            plt.plot(
                x_smooth,
                y_smooth,
                color=colors[idx],
                linestyle=linestyles[idx],
                linewidth=1,  # thinner line,
                label=method_labels[method]
            )
        else:
            plt.plot(
                x_values,
                means,
                '.',  # only markers
                color=colors[idx],
                markersize=0.0000001,
                label=method_labels[method]
            )

    plt.xlabel("Noise Level 10^2")
    plt.ylabel(metric_labels[metric])
    if metric in ["delta_r", "delta_t", "delta_uv", "delta_f"]:
        plt.yscale("log")
        plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, subs=[1.0]))
        plt.gca().yaxis.set_major_formatter(FuncFormatter(customPlotFun))
    if metric in ["delta_skew"]:
        plt.ylim(0, 2)
    plt.title(f"{metric_labels[metric]} vs Noise")
    plt.legend()
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.tight_layout()
    os.makedirs("./synthetic/plots", exist_ok=True)
    plt.savefig(f"./synthetic/plots/{metric}.pdf", dpi=300)

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
