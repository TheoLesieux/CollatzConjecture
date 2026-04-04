"""
Plot the supercycle return mechanism (single panel).

  Ie-to-Ie supercycle return products P = a_{n+L}/a_n by type A/B/C.

Generates: paper/figures/fig_return_mechanism.pdf  (+ .png)
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

from _collatz_common import (
    subband_name, load_data, FIG_DIR,
)

D = load_data()
ns = np.array(D['ns'], dtype=float)
eps_arr = np.array(D['eps'], dtype=float)
an_arr = np.array(D['an'], dtype=float)
N = D['N']

subs_refined = [subband_name(e) for e in eps_arr]

fig, ax = plt.subplots(1, 1, figsize=(8, 5.5))

# ── Supercycle return products ──────────────────────────────
ie_idx = [i for i in range(N) if subs_refined[i] == 'Ie']

sc_A_n, sc_A_P = [], []
sc_B_n, sc_B_P = [], []
sc_C_n, sc_C_P = [], []

for k in range(len(ie_idx) - 1):
    i0, i1 = ie_idx[k], ie_idx[k + 1]
    L  = i1 - i0
    n0 = ns[i0]
    P_obs = an_arr[i1] / an_arr[i0]
    route = tuple(subs_refined[i] for i in range(i0, i1 + 1))
    if L == 7 and 'Id' in route:
        sc_A_n.append(n0); sc_A_P.append(P_obs)
    elif L == 5 and 'Ib' in route:
        sc_B_n.append(n0); sc_B_P.append(P_obs)
    elif L == 5 and 'Ic' in route:
        sc_C_n.append(n0); sc_C_P.append(P_obs)

ax.scatter(sc_A_n, sc_A_P, s=3, color='#2166ac', alpha=0.5,
           label=fr'$\mathcal{{A}}$ ($L=7$, $N={len(sc_A_n)}$)')
ax.scatter(sc_B_n, sc_B_P, s=3, color='#e31a1c', alpha=0.5,
           label=fr'$\mathcal{{B}}$ ($L=5$, $N={len(sc_B_n)}$)')
ax.scatter(sc_C_n, sc_C_P, s=3, color='#33a02c', alpha=0.5,
           label=fr'$\mathcal{{C}}$ ($L=5$, $N={len(sc_C_n)}$)')
ax.axhline(1.0, color='black', linewidth=1.0, zorder=0)

for arr_n, arr_P, color, name in [
    (sc_A_n, sc_A_P, '#2166ac', r'\mathcal{A}'),
    (sc_B_n, sc_B_P, '#e31a1c', r'\mathcal{B}'),
    (sc_C_n, sc_C_P, '#33a02c', r'\mathcal{C}'),
]:
    mask_large = np.array(arr_n) >= 2000
    if mask_large.any():
        mean_P = np.mean(np.array(arr_P)[mask_large])
        ax.axhline(mean_P, color=color, ls='--', lw=0.8,
                   label=f'$P_{{{name}}}^\\infty \\approx {mean_P:.4f}$')

ax.set_xlabel('$n$', fontsize=12)
ax.set_ylabel(r'Supercycle return product $P = a_{n+L}/a_n$', fontsize=12)
ax.set_title(r'Supercycle contraction ($\mathrm{I_e}$-to-$\mathrm{I_e}$)', fontsize=11)
ax.legend(fontsize=8)
ax.set_ylim(0.3, 1.1)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig_return_mechanism.pdf'), dpi=250, bbox_inches='tight')
plt.savefig(os.path.join(FIG_DIR, 'fig_return_mechanism.png'), dpi=250, bbox_inches='tight')
plt.close()
print("Saved fig_return_mechanism.pdf/png")
