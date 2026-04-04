"""
Bilinear regression of log2(a_n) per 3-band and supercycle return products.

Fits:  log2(a_n) = d + c1*n + m*eps(n) + p*log2(n)   per band I/II/III.
Then derives theoretical return products P_A, P_B and verifies against data.
"""
import numpy as np
import mpmath

from _collatz_common import (
    alpha, lam0, band3, load_data,
)

D = load_data()
# numpy arrays only for linalg fitting
ns = np.array(D['ns'], dtype=float)
eps_arr = np.array([float(e) for e in D['eps']])
an_arr = np.array([float(a) for a in D['an']])
log2_an = np.array([float(l) for l in D['log2_an']])
band3_arr = np.array([band3(float(e)) for e in D['eps']])

# ── Global fit for m (shared across all bands) ─────────────
print("=" * 60)
print("BILINEAR FIT: log2(a_n) = d + n*log2(lam0) + m*eps + p*log2(n)")
print("=" * 60)

global_mask = ns >= 200
n_all  = ns[global_mask].astype(float)
e_all  = eps_arr[global_mask]
y_all  = log2_an[global_mask]
X_all  = np.column_stack([np.ones(global_mask.sum()), n_all, e_all, np.log2(n_all)])
coeffs_all, _, _, _ = np.linalg.lstsq(X_all, y_all, rcond=None)
m = coeffs_all[2]
print(f"  Global m = {m:.6f}  (fitted on all bands jointly)\n")

# ── Per-band fit for d, c1, p  (m fixed) ───────────────────
models = {}
for b_id, b_name in enumerate(['I', 'II', 'III']):
    mask = (band3_arr == b_id) & (ns >= 200)
    n_v = ns[mask].astype(float)
    e_v = eps_arr[mask]
    y_v = log2_an[mask] - m * e_v          # absorb m*eps into rhs
    X = np.column_stack([np.ones(mask.sum()), n_v, np.log2(n_v)])
    coeffs, _, _, _ = np.linalg.lstsq(X, y_v, rcond=None)
    pred = X @ coeffs + m * e_v
    resid = log2_an[mask] - pred
    models[b_name] = {'d': coeffs[0], 'c1': coeffs[1], 'p': coeffs[2]}
    print(f"  Band {b_name}: d={coeffs[0]:.4f}, c1={coeffs[1]:.8f} "
          f"(log2(lam0)={float(mpmath.log(lam0)/mpmath.log(2)):.8f}), "
          f"p={coeffs[2]:.4f}, RMSE={np.sqrt(np.mean(resid**2)):.5f}")

# ── Theoretical return products ─────────────────────────────
delta_B = float(2 * alpha - 1)
delta_A = float(3 * alpha - 1)
log2_lam0 = float(mpmath.log(lam0) / mpmath.log(2))
log2_P_B = 2 * log2_lam0 + m * delta_B
log2_P_A = 3 * log2_lam0 + m * delta_A
P_B = 2 ** log2_P_B
P_A = 2 ** log2_P_A

print(f"\nReturn products (asymptotic):")
print(f"  B-cycle (I->III->I): P_inf = {P_B:.6f}, margin {(1-P_B)*100:.1f}%")
print(f"  A-cycle (I->II->III->I): P_inf = {P_A:.6f}, margin {(1-P_A)*100:.1f}%")

# ── Empirical verification ──────────────────────────────────
bandI_idx = [i for i in range(D['N']) if band3_arr[i] == 0 and ns[i] >= 200]
B_prods, A_prods = [], []
for j in range(1, len(bandI_idx)):
    i1, i2 = bandI_idx[j - 1], bandI_idx[j]
    gap = int(ns[i2] - ns[i1])
    P_obs = an_arr[i2] / an_arr[i1]
    if gap == 2:
        B_prods.append(P_obs)
    elif gap == 3:
        A_prods.append(P_obs)

print(f"\nEmpirical check (n >= 200):")
print(f"  B-cycle: count={len(B_prods)}, max={max(B_prods):.6f}, all<1: {all(p<1 for p in B_prods)}")
print(f"  A-cycle: count={len(A_prods)}, max={max(A_prods):.6f}, all<1: {all(p<1 for p in A_prods)}")
