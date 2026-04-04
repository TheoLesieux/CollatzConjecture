"""
Hybrid estimate table (Section 3.4).

Computes the table comparing exact partial sums with the polylogarithm
tail estimate at each convergent denominator q_k of gamma = log2(3).

Columns: k, q_k, exact sum_{n<q_k} a_n, hybrid estimate, digits of 1/2.
"""
import mpmath

from _collatz_common import gamma, alpha, lam0, load_data

C_0 = mpmath.sqrt(2 * mpmath.pi * alpha * gamma)

# Load data
D = load_data()
ns  = D['ns']
eps = D['eps']
an  = D['an']
N   = D['N']

# Compute phi(epsilon) = (gamma/(2*alpha))^epsilon
phi_vals = [(gamma / (2 * alpha)) ** e for e in eps]

# Compute R(n) = a_n * C_0 * n^{3/2} / (lam0^n * phi(epsilon(n)))
R_vals = [an[i] * C_0 * mpmath.mpf(ns[i])**mpmath.mpf('1.5') / (lam0**ns[i] * phi_vals[i])
          for i in range(N)]

# <phi R^0> ~ sample mean of phi(eps(n)) * R(n) for n >= 1000
products = [phi_vals[i] * R_vals[i] for i in range(N) if ns[i] >= 1000]
phi_R0_mean = mpmath.fsum(products) / len(products)
print(f"<phi R^0> = {mpmath.nstr(phi_R0_mean, 10)}  (sample mean for n >= 1000)")
print(f"C_0 = {mpmath.nstr(C_0, 15)}")
print(f"lambda_0 = {mpmath.nstr(lam0, 20)}")
print()

# Polylogarithm tail: Li^{>=N0}_{3/2}(lam0) = sum_{n>=N0} lam0^n / n^{3/2}
def polylog_tail(N0, z, s=mpmath.mpf('1.5'), nmax=50000):
    """Compute sum_{n>=N0} z^n / n^s."""
    total = mpmath.mpf(0)
    for n in range(N0, nmax + 1):
        term = z**n / mpmath.mpf(n)**s
        total += term
        if abs(term) < mpmath.mpf('1e-40'):
            break
    return total

# Total exact sum
total_sum = mpmath.fsum(an)
print(f"Total sum (N=5000) = {mpmath.nstr(total_sum, 30)}")
print()

# Convergent denominators of gamma
convergents = [
    (0,  1),
    (2,  2),
    (3,  5),
    (4,  12),
    (5,  41),
    (6,  53),
    (7,  306),
    (8,  665),
]

print(f"{'k':>3}  {'q_k':>5}  {'exact sum':>20}  {'hybrid':>20}  {'digits':>8}")
print("-" * 65)

for k, qk in convergents:
    # Exact partial sum: sum_{n=1}^{q_k - 1} a_n
    if qk <= 1:
        exact_partial = mpmath.mpf(0)
    else:
        exact_partial = mpmath.fsum(an[:qk - 1])

    # Mean tail: <phi R^0>/C_0 * Li^{>=q_k}_{3/2}(lam0)
    li_tail = polylog_tail(max(qk, 1), lam0)
    mean_tail = phi_R0_mean / C_0 * li_tail

    # Hybrid estimate
    hybrid = exact_partial + mean_tail

    # Digits of 1/2: -log10(|hybrid - 0.5| / 0.5)
    err = abs(hybrid - mpmath.mpf('0.5'))
    if err > 0:
        digits = float(-mpmath.log10(2 * err))
    else:
        digits = float('inf')

    print(f"{k:>3}  {qk:>5}  {mpmath.nstr(exact_partial, 15):>20}  {mpmath.nstr(hybrid, 15):>20}  {digits:>8.1f}")
