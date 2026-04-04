"""
Constant-kernel estimate of u_infty = sum a_n.

Computes:
  - c_+(d), c_-(d)  (carry-split kernels)
  - K(d) = (1-eps(d))*c_+(d) + eps(d)*c_-(d)  (phase-averaged kernel)
  - B = sum b_n  (backbone sum)
  - u_infty ~ B / (1 + sum K(d))  (constant-kernel approximation)

Corresponds to Section 3.3 of the paper (eq. 34).
"""
import mpmath

from _collatz_common import (
    gamma, alpha, beta, lam0, s, eps, load_data
)

DMAX = 500

# ── Compute c_+(d), c_-(d), K(d) ───────────────────────────
c_plus  = [mpmath.mpf(0)] * (DMAX + 1)
c_minus = [mpmath.mpf(0)] * (DMAX + 1)
K_d     = [mpmath.mpf(0)] * (DMAX + 1)

for d in range(1, DMAX + 1):
    sd = s(d)
    ed = eps(d)
    c_plus[d] = mpmath.binomial(sd, d) / mpmath.power(2, sd)
    if sd - 1 >= d:
        c_minus[d] = mpmath.binomial(sd - 1, d) / mpmath.power(2, sd - 1)
    K_d[d] = (1 - ed) * c_plus[d] + ed * c_minus[d]

P_plus = mpmath.fsum(c_plus[1:])
S_K    = mpmath.fsum(K_d[1:])

# ── Backbone sum B = sum b_n ────────────────────────────────
D = load_data()
B_sum = mpmath.fsum(D['bn'])

# ── Exact partial sum ───────────────────────────────────────
S_exact = mpmath.fsum(D['an'])

# ── Constant-kernel estimate ────────────────────────────────
# Summing the renewal equation: S = B - (sum K) * S ⟹ S = B/(1+sum K)
S_Kavg = B_sum / (1 + S_K)

# ── Report ──────────────────────────────────────────────────
print("=" * 60)
print("CONSTANT-KERNEL ESTIMATE  (Section 3.3)")
print("=" * 60)
print(f"gamma      = {float(gamma):.12f}")
print(f"alpha      = {float(alpha):.12f}")
print(f"beta       = {float(beta):.12f}")
print(f"lambda_0   = {float(lam0):.12f}")
print()
print(f"P_+ = sum c_+(d)            = {mpmath.nstr(P_plus, 15)}")
print(f"1 + P_+                     = {mpmath.nstr(1 + P_plus, 15)}")
print(f"sum K(d)  [phase-averaged]  = {mpmath.nstr(S_K, 15)}")
print(f"1 + sum K                   = {mpmath.nstr(1 + S_K, 15)}")
print(f"B = sum b_n                 = {mpmath.nstr(B_sum, 15)}")
print()
print(f"S_exact = sum a_n (N=5000)  = {mpmath.nstr(S_exact, 30)}")
print(f"u_infty ~ B/(1+sum K)       = {mpmath.nstr(S_Kavg, 15)}")
print(f"  error from 1/2            = {mpmath.nstr(abs(S_Kavg - mpmath.mpf('0.5')), 6)}  ({mpmath.nstr(abs(S_Kavg - mpmath.mpf('0.5'))/mpmath.mpf('0.5')*100, 4)}%)")
