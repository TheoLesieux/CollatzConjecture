"""
Plot the 11-subband structure of rho(n) = a_{n+1}/a_n vs Sturmian phase.

Generates: paper/figures/fig_rho_subbands.pdf  (+ .png)
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from _collatz_common import (
    gamma, alpha, beta, jbeta,
    thresholds_11, sub_names, sub_colors,
    load_data, FIG_DIR,
)
import os

D = load_data()
ns = np.array(D['ns'], dtype=float)
eps_arr = np.array(D['eps'], dtype=float)
an_arr = np.array(D['an'], dtype=float)
rho   = an_arr[1:] / an_arr[:-1]

fig, ax = plt.subplots(figsize=(10, 5.5))
nmin = 500
mask_n = ns[:-1] >= nmin

thr_info = [
    (jbeta(5),  r'$\{5\beta\}$'),  (jbeta(10), r'$\{10\beta\}$'),
    (2*gamma-3, r'$2\gamma{-}3$'), (jbeta(3),  r'$\{3\beta\}$'),
    (beta,      r'$\beta$'),       (alpha,     r'$\alpha$'),
    (jbeta(4),  r'$\{4\beta\}$'),  (jbeta(2),  r'$2\beta$'),
    (jbeta(7),  r'$\{7\beta\}$'),  (jbeta(12), r'$\{12\beta\}$'),
]

for i in range(len(thresholds_11) - 1):
    lo, hi = thresholds_11[i], thresholds_11[i + 1]
    mask = mask_n & (eps_arr[:-1] >= lo) & (eps_arr[:-1] < hi)
    if mask.sum() == 0:
        continue
    e = eps_arr[:-1][mask]
    r = rho[mask]
    ax.scatter(e, r, s=2, alpha=0.25, color=sub_colors[i], rasterized=True, zorder=2)
    rmean = r.mean()
    ax.hlines(rmean, lo + 0.002, hi - 0.002, colors=sub_colors[i], linewidths=3.0, zorder=6)
    emid = (lo + hi) / 2
    ax.text(emid, rmean + 0.015, sub_names[i], ha='center', va='bottom',
            fontsize=8, fontweight='bold', color=sub_colors[i], zorder=7)

for t, label in thr_info:
    is_main = abs(t - beta) < 0.01 or abs(t - alpha) < 0.01
    ax.axvline(t, color='0.3' if is_main else '0.6',
               linewidth=1.2 if is_main else 0.6,
               linestyle='-' if is_main else ':', zorder=1)

ax.set_xlabel(r'Phase $\varepsilon(n)$', fontsize=12)
ax.set_ylabel(r'$\rho(n) = a_{n+1}/a_n$', fontsize=12)
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(0.55, 1.72)
ax.set_title(r'Sub-band structure of $\rho(\varepsilon)$', fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig_rho_subbands.pdf'), dpi=250, bbox_inches='tight')
plt.savefig(os.path.join(FIG_DIR, 'fig_rho_subbands.png'), dpi=250, bbox_inches='tight')
plt.close()
print("Saved fig_rho_subbands.pdf/png")
