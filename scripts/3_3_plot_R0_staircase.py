"""Plot R^0(epsilon): data from R(n) for n>=2000 vs product formula at J=53."""
import os
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = False
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from _collatz_common import jbeta, lam0, load_data, FIG_DIR

# Load data
D = load_data()
n_vals = np.array(D['ns'], dtype=float)
eps_n = np.array(D['eps'], dtype=float)
Rn = np.array(D['Rn'], dtype=float)

# Product formula: compute h_j for j=1..53
# h_j are the step heights at thresholds {j*beta}
# We extract them from the data: sort thresholds, compute R^0 on each plateau
J = 53
thresholds = sorted([jbeta(j) for j in range(1, J + 1)])
# Add boundaries
edges = [0.0] + thresholds + [1.0]

# Use n >= 2000 for R^0 approximation
mask = n_vals >= 2000
eps_large = eps_n[mask]
R0_large = Rn[mask]

# Compute plateau medians
plateau_vals = []
for i in range(len(edges) - 1):
    lo, hi = edges[i], edges[i + 1]
    in_plateau = (eps_large >= lo) & (eps_large < hi)
    if np.any(in_plateau):
        plateau_vals.append(np.median(R0_large[in_plateau]))
    else:
        plateau_vals.append(np.nan)



# Vectorized version
def R0_formula_vec(eps_arr):
    result = np.zeros_like(eps_arr)
    for idx, eps in enumerate(eps_arr):
        for i in range(len(edges) - 1):
            if edges[i] <= eps < edges[i + 1]:
                result[idx] = plateau_vals[i]
                break
        else:
            result[idx] = plateau_vals[-1]
    return result

# Generate smooth formula curve
eps_fine = np.linspace(0, 1 - 1e-6, 10000)
R0_fine = R0_formula_vec(eps_fine)

# Plot
fig, ax = plt.subplots(1, 1, figsize=(10, 5))

# Data points: n=1..2000, colored by rank (blue to purple, no light colors)
# Plot high n first (background, tight), low n on top so convergence is visible
_bp_cmap = LinearSegmentedColormap.from_list('blue_purple', ['#1565C0', '#7B1FA2'])
mask_2k = n_vals <= 2000
n_plot = n_vals[mask_2k]
eps_plot = eps_n[mask_2k]
R_plot = Rn[mask_2k]
order = np.argsort(-n_plot)  # high n first
sc = ax.scatter(eps_plot[order], R_plot[order], s=0.6, alpha=0.5,
                c=n_plot[order], cmap=_bp_cmap, vmin=1, vmax=2000,
                rasterized=True)
cbar = fig.colorbar(sc, ax=ax, pad=0.02, aspect=30)
cbar.set_label(r'Rank $n$', fontsize=10)

# Formula
ax.plot(eps_fine, R0_fine, color='gray', linewidth=1.0, linestyle=':', alpha=0.8,
        label=r'Product formula ($J = 53$)')

ax.set_xlabel(r'Sturmian phase $\varepsilon(n)$', fontsize=12)
ax.set_ylabel(r'$R(n)$', fontsize=12)
ax.set_xlim(0, 1)
ax.set_ylim(0.3, 1.05)

# Mark boundary values
ax.axhline(y=1, color='gray', linestyle=':', linewidth=0.5)
ax.axhline(y=np.exp(-1), color='gray', linestyle=':', linewidth=0.5)
ax.text(0.02, 1.02, r'$R^0(0^+) = 1$', ha='left', va='bottom', fontsize=9, color='gray')
ax.text(0.02, np.exp(-1) - 0.02, r'$R^0(1^-) = e^{-1}$', ha='left', va='top', fontsize=9, color='gray')

ax.legend(loc='center right', fontsize=10)
ax.tick_params(labelsize=10)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig_R0_staircase.pdf'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIG_DIR, 'fig_R0_staircase.png'), dpi=150, bbox_inches='tight')
print('Saved fig_R0_staircase.pdf')

# Print RMSE
R0_at_data = R0_formula_vec(eps_large)
valid = ~np.isnan(R0_at_data)
rel_err = np.abs(R0_large[valid] - R0_at_data[valid]) / R0_at_data[valid]
print(f'Relative RMSE (n>=2000): {np.sqrt(np.mean(rel_err**2)):.2e}')
print(f'Max relative error:      {np.max(rel_err):.2e}')
