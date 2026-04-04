"""Plot Phi(epsilon(D)) in batches of q6=53, showing how each 'period'
traces the same eigenfunction curve.

Uses the ALREADY-COMPUTED data from inductive_bootstrap.py approach:
  We recompute V^[D] incrementally to avoid O(D_MAX * N * D) cost.
  Key trick: V^[D] only differs from V^[D-1] for n >= D+1,
  and to get excess(D) we only need sum(V^[D][n] * 2^{S-s_n}).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import comb
import mpmath, os, time

from _collatz_common import gamma, s, lam0, FIG_DIR

N = 1200
D_MAX = 500

CACHE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "excess_D0_500.txt")

if os.path.exists(CACHE):
    print("Loading cached excess data...")
    excess = {}
    with open(CACHE) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            d_str, e_str = line.split()
            excess[int(d_str)] = float(e_str)
else:
    print("Computing excess(D) for D=0..500 ...", flush=True)
    t0 = time.time()
    sn = [s(n) for n in range(N + 1)]
    S = sn[N]
    backbone = [0] * (N + 1)
    for n in range(1, N + 1):
        backbone[n] = comb(sn[n] - 1, n - 1)
    half_num = 1 << (S - 1)

    u0_int = sum(backbone[n] * (1 << (S - sn[n])) for n in range(1, N + 1))
    excess_int = {0: u0_int - half_num}

    for D in range(1, D_MAX + 1):
        V = [0] * (N + 1)
        for n in range(1, N + 1):
            total = backbone[n]
            for d2 in range(1, min(D, n - 1) + 1):
                diff = sn[n] - sn[n - d2]
                b = comb(diff, d2) if diff >= d2 else 0
                total -= b * V[n - d2]
            V[n] = total
        u_int = sum(V[n] * (1 << (S - sn[n])) for n in range(1, N + 1))
        excess_int[D] = u_int - half_num
        if D % 50 == 0:
            print(f"  D={D}  ({time.time()-t0:.0f}s)", flush=True)

    S_val = S
    excess = {}
    for D in range(0, D_MAX + 1):
        excess[D] = float(mpmath.mpf(excess_int[D]) / mpmath.power(2, S_val))

    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    with open(CACHE, 'w') as f:
        f.write("# D  excess(D) = (u^[D] - 1/2)  for N=1200, D=0..500\n")
        f.write("# Columns: D  excess(D)\n")
        for D_key in range(0, D_MAX + 1):
            f.write(f"{D_key}\t{excess[D_key]:.17e}\n")
    print(f"Cached to {CACHE}  ({time.time()-t0:.0f}s total)")

# Build arrays
D_all = np.arange(1, D_MAX + 1)
eps_all = np.array([float(mpmath.ceil(gamma * D) - gamma * D) for D in D_all])
phi_all = np.array([excess[D] / float(lam0**D) for D in D_all])

# ── Figure: color gradient by rank D (like R0 staircase) ──
fig, ax = plt.subplots(figsize=(10, 6))

from matplotlib.colors import LinearSegmentedColormap
_bp_cmap = LinearSegmentedColormap.from_list('blue_purple', ['#1565C0', '#7B1FA2'])

# Plot high D first (background), low D on top so convergence is visible
order = np.argsort(-D_all)
sc = ax.scatter(eps_all[order], phi_all[order], s=4, alpha=0.6,
                c=D_all[order], cmap=_bp_cmap, vmin=1, vmax=D_MAX,
                rasterized=True)
cbar = fig.colorbar(sc, ax=ax, pad=0.02, aspect=30)
cbar.set_label(r'Depth $D$', fontsize=10)

ax.set_xlabel(r'$\varepsilon(D) = \lceil\gamma D\rceil - \gamma D$', fontsize=13)
ax.set_ylabel(r'$\Phi(D) = \mathrm{excess}(D)\,/\,\lambda_0^D$', fontsize=13)
ax.set_title(r'Eigenfunction $\Phi(\varepsilon)$', fontsize=14)
ax.axhline(0, color='red', ls='--', lw=0.8)
ax.set_xlim(-0.02, 1.02)

plt.tight_layout()
os.makedirs(FIG_DIR, exist_ok=True)
out_png = os.path.join(FIG_DIR, 'phi_eigenfunction.png')
out_pdf = os.path.join(FIG_DIR, 'phi_eigenfunction.pdf')
plt.savefig(out_png, dpi=200, bbox_inches='tight')
plt.savefig(out_pdf, bbox_inches='tight')
print(f"Saved {out_png} and {out_pdf}")
plt.close()
