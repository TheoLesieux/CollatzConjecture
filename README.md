# Collatz Conjecture — Sturmian Analysis of First-Decrease Densities

This repository accompanies the paper *"Exact first-decrease densities for the Collatz map and their Sturmian fine structure"* (Th. Lesieux), which develops a combinatorial framework for the Collatz conjecture based on residue-class pair exploration.

## Repository Structure

```
residue_tree/     Residue-class tree explorer (Python BFS)
scripts/          Figure generation and computational verification
data/             Precomputed sequence data
```

### `residue_tree/`

An independent constructive exploration of the Collatz map at the level of residue classes. Starting from the two root classes (2, 0) and (2, 1), a queue-based BFS applies the reduced Collatz map to pairs (M, r), splitting a class into two children whenever parity depends on the quotient k. Each leaf proves first decrease for an entire infinite arithmetic progression. The algorithm never invokes the analytic recurrence, providing both an independent verification of the densities aₙ and the only available analysis of the k = 0 outcast condition.

```bash
cd residue_tree
python main.py --steps 200000
```

### `scripts/`

Shared helpers live in `_collatz_common.py` (constants γ, α, β, λ₀; data loading; band classifiers). All numerical computation uses `mpmath` at 50-digit precision. Every other script has a single purpose:

| Script | Purpose |
|---|---|
| `_generate_data.py` | Generate `data/collatz_data.txt` |
| `2_5_remainder_decrease_verification.py` | Verify (3/2)^p − 1 < 2^n − 3^p and 2^ε ≥ 1 + 2^{−p} for p = 1…5000 (Thm 2.5) |
| `3_3_constant_kernel_estimate.py` | Compute c±(d), K(d), B = ∑bₙ and u∞ ≈ B/(1+∑K) (§3.3) |
| `3_3_plot_R0_staircase.py` | Plot R⁰(ε) staircase vs product formula at J = 53 |
| `3_4_hybrid_estimate_table.py` | Polylogarithm tail hybrid estimate at CF convergents (§3.4) |
| `3_5_bilinear_fit.py` | Bilinear regression log₂(aₙ) per 3-band + theoretical return products P_A, P_B |
| `3_5_plot_rho_subbands.py` | Figure: ρ(n) = aₙ₊₁/aₙ vs ε(n) — 11-subband structure |
| `3_5_plot_return_mechanism.py` | Figure: supercycle return mechanism |
| `3_5_supercycle_statistics.py` | Supercycle classification (A/B/C), word structure, tail bracket |
| `3_6_exact_verification.py` | Exact integer proof that u^[D] > ½ for D = 1…200 |
| `3_6_kernel_truncation.py` | 50-digit mpmath verification of u^[D] > ½ for 18 selected D |
| `3_6_plot_phi_eigenfunction.py` | Plot eigenfunction Φ(ε) = excess(D)/λ₀^D, colored by depth D (D = 1…500) |
| `A_k0_threshold_plot.py` | Appendix A: k = 0 failure threshold r₀ vs τ figure |
| `A_route_comparison_table.py` | Appendix A: recurrence vs BFS density comparison (p = 0…18, 20M steps) |
| `B_step_heights_plot.py` | Appendix B: step heights h_j of the Devil's staircase R⁰(ε) |

### `data/`

| File | Description |
|---|---|
| `collatz_data.txt` | 5000 terms — columns: n, ε(n), s(n), aₙ, log₂(aₙ), R(n), bₙ, vₙ, Bₙ (50-digit precision) |
| `excess_D0_500.txt` | Cache: excess(D) = u^[D] − ½ for D = 0…500. Regenerated automatically by `3_6_plot_phi_eigenfunction.py` if missing. |

## Requirements

- Python 3.10+
- `mpmath`, `matplotlib`, `numpy`

## License

MIT — see [LICENSE](LICENSE).
