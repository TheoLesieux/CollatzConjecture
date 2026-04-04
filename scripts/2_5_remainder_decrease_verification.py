"""
Verify the remainder-decrease inequality (Theorem 2.5) for p = 1..5000.

Three checks per rank p (with n = ceil(gamma * p)):

  (A) Original inequality (eq:rem-ineq):
        (3/2)^p - 1  <  2^n - 3^p

  (B) Sufficient condition (eq:suff-cond):
        2^{epsilon(p)}  >=  1 + 2^{-p}
      where epsilon(p) = ceil(gamma*p) - gamma*p

  (C) Convexity bound:
        epsilon(p) * ln(2)  >=  2^{-p}

All arithmetic uses mpmath at 50-digit precision for exact verification.

Usage:
    python scripts/2_5_remainder_decrease_verification.py
"""
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

import mpmath
from _collatz_common import gamma

mpmath.mp.dps = 50
ln2 = mpmath.log(2)

P_MAX = 5000


def main():
    print(f"Verifying remainder-decrease inequalities for p = 1..{P_MAX}")
    print(f"  Precision: {mpmath.mp.dps} digits\n")

    all_A = True
    all_B = True
    all_C = True
    worst_margin_A = None
    worst_margin_C = None

    for p in range(1, P_MAX + 1):
        n = int(mpmath.ceil(gamma * p))
        eps = n - gamma * p                  # epsilon(p) > 0

        # (A) (3/2)^p - 1  <  2^n - 3^p
        lhs_A = (mpmath.mpf(3) / 2) ** p - 1
        rhs_A = mpmath.mpf(2) ** n - mpmath.mpf(3) ** p
        ok_A = lhs_A < rhs_A
        margin_A = float(rhs_A - lhs_A)

        # (B) 2^eps >= 1 + 2^{-p}
        lhs_B = mpmath.power(2, eps)
        rhs_B = 1 + mpmath.power(2, -p)
        ok_B = lhs_B >= rhs_B

        # (C) eps * ln2 >= 2^{-p}
        lhs_C = eps * ln2
        rhs_C = mpmath.power(2, -p)
        ok_C = lhs_C >= rhs_C
        margin_C = float(mpmath.log10(lhs_C / rhs_C)) if rhs_C > 0 else float('inf')

        if not ok_A:
            print(f"  FAIL (A) at p={p}: LHS={lhs_A}, RHS={rhs_A}")
            all_A = False
        if not ok_B:
            print(f"  FAIL (B) at p={p}")
            all_B = False
        if not ok_C:
            print(f"  FAIL (C) at p={p}")
            all_C = False

        if worst_margin_A is None or margin_A < worst_margin_A[1]:
            worst_margin_A = (p, margin_A)
        if worst_margin_C is None or margin_C < worst_margin_C[1]:
            worst_margin_C = (p, margin_C)

        if p % 1000 == 0:
            print(f"  p = {p:5d}  eps = {float(eps):.6f}  "
                  f"margin(A) = {margin_A:.4e}  "
                  f"log10-margin(C) = {margin_C:.1f}")

    print(f"\nResults (p = 1..{P_MAX}):")
    print(f"  (A) (3/2)^p - 1 < 2^n - 3^p :  {'PASS' if all_A else 'FAIL'}")
    print(f"  (B) 2^eps >= 1 + 2^(-p)      :  {'PASS' if all_B else 'FAIL'}")
    print(f"  (C) eps*ln2 >= 2^(-p)         :  {'PASS' if all_C else 'FAIL'}")
    if not all_B or not all_C:
        print("  Note: (B) and (C) are sufficient but not necessary conditions.")
        print("  They may fail at small p (e.g. p=1) where the original")
        print("  inequality (A) nevertheless holds by direct evaluation.")
    p_w, m_w = worst_margin_A
    print(f"\n  Tightest margin (A): p = {p_w}, "
          f"2^n - 3^p - ((3/2)^p - 1) = {m_w:.4e}")
    p_w, m_w = worst_margin_C
    print(f"  Tightest margin (C): p = {p_w}, "
          f"log10(eps*ln2 / 2^{{-p}}) = {m_w:.1f} orders of magnitude")


if __name__ == "__main__":
    main()
