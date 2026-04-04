"""
Plot k=0 failure threshold analysis (Appendix A, Figure).

For each BFS first-descent leaf at rank p, computes the carry term
C(b) = 2^{s(p)} * r' - 3^p * r_0 and the failure threshold
tau = C(b) / (2^{s(p)} - 3^p).  Plots r_0 vs tau for all leaves.

Usage:
    python scripts/A_k0_threshold_plot.py
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'residue_tree'))
from classes.pair_handler import PairHandler

from _collatz_common import s, load_data, FIG_DIR

BFS_STEPS = 5_000_000


def main():
    # --- Run BFS and collect leaf data ---
    handler = PairHandler()
    leaf_data = defaultdict(list)  # s(p) -> [(r0, r_n, M0)]
    original_step = handler.step

    def patched_step():
        if not handler.queue:
            return False
        pair = handler.queue[0]
        if handler._has_decreased(pair):
            a0, b0 = pair.seed
            at, bt = pair.transformed
            n = a0.bit_length() - 1
            leaf_data[n].append((b0, bt, a0))
        return original_step()

    handler.step = patched_step
    print(f"Running BFS ({BFS_STEPS:,} steps)...", file=sys.stderr)
    handler.run(max_steps=BFS_STEPS)

    # --- Find completed ranks ---
    _d = load_data()
    v = [1] + _d['vn'][:19]  # v[0]=1, v[1..19] from data

    completed = []
    for p in range(20):
        sp = s(p)
        bfs_density = handler.densities.get(sp, 0.0)
        rec_density = v[p] / 2 ** sp
        if abs(bfs_density - rec_density) < 1e-15:
            completed.append(p)
        else:
            break

    p_max = completed[-1]
    print(f"Completed ranks: 0..{p_max}", file=sys.stderr)

    # --- Compute thresholds ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left panel: per-rank scatter of r_0 / tau (safety ratio)
    ax1 = axes[0]
    # Use tab20 for maximally distinguishable colors across ranks
    cmap = plt.cm.tab20
    rank_colors = {p: cmap(p / max(p_max, 1)) for p in range(p_max + 1)}

    tightest_margin = float('inf')
    tightest_info = None
    total_leaves = 0

    rank_data = {}  # p -> (rs_0, taus)

    for p in range(1, p_max + 1):
        sp = s(p)
        data = leaf_data.get(sp, [])
        if not data:
            continue

        gap = 2**sp - 3**p
        rs_0 = []
        taus = []

        for b0, bn, a0 in data:
            if b0 == 0:
                continue
            total_leaves += 1
            C = 2**sp * bn - 3**p * b0
            tau = C / gap if gap > 0 else 0
            rs_0.append(b0)
            taus.append(tau)

            if tau > 0:
                ratio = b0 / tau
                if ratio < tightest_margin:
                    tightest_margin = ratio
                    tightest_info = (p, b0, bn, tau, ratio)

        rank_data[p] = (np.array(rs_0), np.array(taus))

    # Left: scatter r_0 vs tau per rank, log-log
    for p in sorted(rank_data):
        b0s, ts = rank_data[p]
        pos = ts > 0
        if pos.any():
            ax1.scatter(ts[pos], b0s[pos], s=2, alpha=0.25,
                        color=rank_colors[p], rasterized=True,
                        label=f'p={p}' if p % 3 == 1 else '')

    lo, hi = ax1.get_xlim()
    span = np.linspace(max(lo, 0.1), hi, 100)
    ax1.plot(span, span, 'r--', lw=1.2, label=r'$r^{(0)} = \tau$ (failure)')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'Failure threshold $\tau = C(\mathbf{b})/(2^{s(p)}-3^p)$')
    ax1.set_ylabel(r'Seed remainder $r^{(0)}$')
    ax1.set_title(r'$r^{(0)}$ vs $\tau$ (all leaves)')
    ax1.legend(fontsize=7, ncol=2, loc='upper left')

    # Right panel: per-rank minimum safety ratio
    ax2 = axes[1]
    ranks = []
    min_ratios = []

    for p in range(1, p_max + 1):
        sp = s(p)
        b0s, ts = rank_data.get(p, (np.array([]), np.array([])))
        if len(b0s) == 0:
            continue

        pos = ts > 0
        if not pos.any():
            continue

        ratios = b0s[pos] / ts[pos]
        ranks.append(p)
        min_ratios.append(ratios.min())

    ax2.bar(ranks, min_ratios, color=[rank_colors[p] for p in ranks], edgecolor='k',
            linewidth=0.3)
    ax2.axhline(1, color='r', ls='--', lw=1, label=r'Failure line $r^{(0)}/\tau = 1$')
    ax2.set_xlabel('Rank $p$')
    ax2.set_ylabel(r'Minimum $r^{(0)}/\tau$ across leaves')
    ax2.set_title('Safety margin by rank')
    ax2.legend(fontsize=9)
    ax2.set_yscale('log')

    fig.tight_layout()
    os.makedirs(FIG_DIR, exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(FIG_DIR, f'fig_k0_threshold.{ext}'), dpi=200)

    print(f"\nTotal leaves (r_0 > 0): {total_leaves}", file=sys.stderr)
    if tightest_info:
        p, b0, bn, tau, ratio = tightest_info
        print(f"Tightest: p={p}, r_0={b0}, r_n={bn}, "
              f"tau={tau:.2f}, ratio={ratio:.4f}", file=sys.stderr)
    print(f"Saved to {FIG_DIR}/fig_k0_threshold.pdf", file=sys.stderr)


if __name__ == "__main__":
    main()
