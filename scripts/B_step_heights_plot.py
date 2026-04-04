"""Plot step heights h_j of the Devil's staircase R^0(eps) = lam0^{H(eps)}."""
import os
import numpy as np
import matplotlib.pyplot as plt

from _collatz_common import jbeta, lam0, load_data, FIG_DIR

# Load data
D = load_data()
n_vals = np.array(D['ns'], dtype=float)
eps_n = np.array(D['eps'], dtype=float)
Rn = np.array(D['Rn'], dtype=float)

J = 53
thresh_j = [jbeta(j) for j in range(1, J + 1)]

# Sort thresholds keeping track of original j
sorted_pairs = sorted(enumerate(thresh_j, 1), key=lambda x: x[1])
sorted_j = [p[0] for p in sorted_pairs]
sorted_thresh = [p[1] for p in sorted_pairs]
edges = [0.0] + sorted_thresh + [1.0]

# Use n >= 2000 for stable plateau medians
mask = n_vals >= 2000
eps_large = eps_n[mask]
R0_large = Rn[mask]

# Compute plateau medians -> H values
plateau_H = []
for i in range(len(edges) - 1):
    lo, hi = edges[i], edges[i + 1]
    in_plateau = (eps_large >= lo) & (eps_large < hi)
    if np.any(in_plateau):
        med = np.median(R0_large[in_plateau])
        plateau_H.append(np.log(med) / np.log(float(lam0)))
    else:
        plateau_H.append(np.nan)

# Step heights h_j (mapped back to original j index)
h_by_j = {}
for i in range(len(sorted_j)):
    j = sorted_j[i]
    h_by_j[j] = plateau_H[i + 1] - plateau_H[i]

js = np.arange(1, J + 1, dtype=float)
hs = np.array([h_by_j[j] for j in range(1, J + 1)])

# Power-law fit for j >= 2
mask_fit = js >= 2
log_j = np.log(js[mask_fit])
log_h = np.log(hs[mask_fit])
theta, log_zeta = np.polyfit(log_j, log_h, 1)
zeta = np.exp(log_zeta)
print(f'Power-law fit: h_j ~ {zeta:.4f} * j^({theta:.4f})')
print(f'Sum h_j = {np.sum(hs):.4f}')

# Plot
fig, ax = plt.subplots(figsize=(8, 4.5))

ax.semilogy(js, hs, 'o', color='#1565C0', markersize=4.5, zorder=3,
            label=r'$h_j$ (data)')

# Power-law fit line
j_fit = np.linspace(1, 53, 200)
ax.semilogy(j_fit, zeta * j_fit**theta, '-', color='#7B1FA2', linewidth=1.5,
            label=rf'$h_j \approx {zeta:.2f}\,j^{{{theta:.2f}}}$')

ax.set_xlabel(r'$j$', fontsize=12)
ax.set_ylabel(r'$h_j$', fontsize=12)
ax.set_xlim(0, 55)
ax.legend(fontsize=11, framealpha=0.9)
ax.tick_params(labelsize=10)
ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig_step_heights.pdf'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(FIG_DIR, 'fig_step_heights.png'), dpi=150, bbox_inches='tight')
print('Saved fig_step_heights.pdf')
