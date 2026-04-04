"""
Generate the route-comparison table (Appendix A, Table 1).

Computes a_p via two independent routes:
  (1) Prime combination recurrence (Corollary 2.8)
  (2) Breadth-first residue-class tree exploration
and verifies they agree at every completed rank.

Usage:
    python scripts/A_route_comparison_table.py [--steps 20000000]
"""
import argparse
import sys
import os

# ── Import BFS engine from residue_tree ─────────────────────
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'residue_tree'))
from classes.pair_handler import PairHandler

from _collatz_common import s, load_data

P_MAX = 18  # upper bound; table stops at last completed rank


def main():
    parser = argparse.ArgumentParser(description="Route-comparison table (Appendix A)")
    parser.add_argument("--steps", type=int, default=20_000_000,
                        help="Maximum BFS iterations (default: 20000000)")
    args = parser.parse_args()

    # Route 1: from precomputed data
    _d = load_data()
    v = [1] + _d['vn'][:P_MAX]  # v[0]=1, v[1..P_MAX] from data
    a_rec = {p: v[p] / 2 ** s(p) for p in range(P_MAX + 1)}

    # Route 2: BFS
    print(f"Running BFS ({args.steps:,} steps)...", file=sys.stderr)
    handler = PairHandler()
    handler.run(max_steps=args.steps)

    # Map BFS densities (keyed by step count n) to ranks p
    a_bfs = {}
    for p in range(P_MAX + 1):
        a_bfs[p] = handler.densities.get(s(p), 0.0)

    # Find last fully completed rank
    last_complete = -1
    for p in range(P_MAX + 1):
        if abs(a_rec[p] - a_bfs[p]) < 1e-15:
            last_complete = p
        else:
            break

    # Print table
    cumulative = 0.0
    print(f"\n{'p':>3}  {'s(p)':>4}  {'a_p (recurrence)':>22}  "
          f"{'a_p (BFS)':>22}  {'Match':>5}")
    print("-" * 64)
    for p in range(last_complete + 2):  # show one incomplete rank
        if p > P_MAX:
            break
        cumulative += a_rec[p]
        match = abs(a_rec[p] - a_bfs[p]) < 1e-15
        mark = 'Y' if match else '.'
        print(f"{p:>3}  {s(p):>4}  {a_rec[p]:>22.15f}  "
              f"{a_bfs[p]:>22.15f}  {mark}")
    print("-" * 64)
    cum_complete = sum(a_rec[p] for p in range(last_complete + 1))
    print(f"{'':>3}  {'':>4}  {cum_complete:>22.15f}")
    print(f"\nCompleted ranks: 0..{last_complete}  "
          f"({args.steps:,} BFS steps)")


if __name__ == "__main__":
    main()
