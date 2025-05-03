import matplotlib.pyplot as plt
import numpy as np

# Static data
error_types = [r'$\Delta f$', r'$\Delta c$']
methods = ['Ours', 'Right Cylinder', 'Gummeson']

errors = np.array([
    [0.06, 0.32],   # Ours: Δf, Δc
    [np.nan, np.nan],  # Right Cylinder: no valid solution
    [0.0, np.nan]    # Gummeson: Δf, Δc
])

# Parameters
x = np.arange(len(error_types))  # X locations for error types
width = 0.25  # Width of the bars
nan_display_value = 1  # Fixed value for NaNs

# Bar colors for each method
colors = ['skyblue', 'lightgray', 'salmon']

fig, ax = plt.subplots()

# Prepare bars
for method_idx, method in enumerate(methods):
    offset = (method_idx - 1) * width
    for error_idx in range(len(error_types)):
        error_value = errors[method_idx, error_idx]
        if np.isnan(error_value):
            # Show a hatched bar for NaN
            rect = ax.bar(
                x[error_idx] + offset, nan_display_value, width,
                color=colors[method_idx],
                hatch='//',
                edgecolor='black'
            )
            ax.annotate('N/A',
                        xy=(x[error_idx] + offset, nan_display_value),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom')
        else:
            # Normal bar
            rect = ax.bar(
                x[error_idx] + offset, error_value, width,
                color=colors[method_idx]
            )
            ax.annotate(f'{error_value:.1f}',
                        xy=(x[error_idx] + offset, error_value),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom')

# Customize plot
ax.set_ylabel('Error')
ax.set_title('Estimation Errors on Real Data')
ax.set_xticks(x)
ax.set_xticklabels(error_types)
ax.legend(methods)

# Add grid
ax.yaxis.grid(True, linestyle='--', which='major', color='grey', alpha=0.5)

fig.tight_layout()
plt.savefig('error_reals.pdf')
plt.show()